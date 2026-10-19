"""Deterministic, generated BAM/FRAG benchmark inputs; no external downloads."""

import hashlib
import json
from pathlib import Path
import random
import struct
import zlib


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class BGZFWriter:
    def __init__(self, path):
        self.handle = Path(path).open('wb')
        self.pending = bytearray()

    def block(self, data):
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = compressor.compress(data) + compressor.flush()
        size = len(payload) + 26
        assert size <= 65536
        self.handle.write(b'\x1f\x8b\x08\x04' + b'\0' * 4 + b'\0\xff'
                          + struct.pack('<H2sHH', 6, b'BC', 2, size - 1)
                          + payload + struct.pack('<II', zlib.crc32(data), len(data)))

    def write(self, data):
        self.pending.extend(data)
        while len(self.pending) >= 60000:
            self.block(bytes(self.pending[:60000]))
            del self.pending[:60000]

    def close(self):
        if self.pending:
            self.block(bytes(self.pending))
        self.block(b'')
        self.handle.close()


def bam_record(chrom, left, index):
    name = f'r{index}'.encode() + b'\0'
    # Complete pair; both reads contribute to BAM, one interval to BAMPE.
    for flag, pos, mate, length in ((99, left, left + 114, 150),
                                    (147, left + 114, left, -150)):
        bin_id = next((offset + (pos >> shift) for shift, offset in
                       ((14, 4681), (17, 585), (20, 73), (23, 9), (26, 1))
                       if pos >> shift == (pos + 35) >> shift), 0)
        body = (struct.pack('<iiIIiiii', chrom, pos,
                            (bin_id << 16) | (30 << 8) | len(name),
                            (flag << 16) | 1, 36, chrom, mate, length)
                + name + struct.pack('<I', 36 << 4)
                + b'\x12' * 18 + b'\x1e' * 36)
        yield struct.pack('<i', len(body)) + body


def prepare(directory, small=6000, large=500000, contigs=256):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if min(small, large, contigs) < 1:
        raise ValueError('record counts and contigs must be positive')
    cases = []
    for size, count, n_chrom in [('small', small, 3), ('large', large, 3),
                                 ('many', large, contigs)]:
        inputs = {'BAM': {}, 'FRAG': {}}
        for role in ('treat', 'control'):
            paths = {kind: directory / f'{size}-{role}.{suffix}'
                     for kind, suffix in [('BAM', 'bam'), ('FRAG', 'tsv.gz')]}
            bam, frag = BGZFWriter(paths['BAM']), BGZFWriter(paths['FRAG'])
            header = b'BAM\1' + struct.pack('<ii', 0, n_chrom)
            for chrom in range(n_chrom):
                name = f'chr{chrom:04}'.encode() + b'\0'
                header += struct.pack('<i', len(name)) + name + struct.pack('<i', 2000000)
            bam.write(header)
            rng = random.Random(773 + (role == 'control'))
            index = 0
            for chrom in range(n_chrom):
                number = count // n_chrom + (chrom < count % n_chrom)
                # Several peaks plus a dispersed background; unique positions
                # and barcode keys exercise more than duplicate compression.
                positions = sorted(
                    (100000 + (i % 16) * 70000 + rng.randrange(600)
                     if role == 'treat' and i % 4 else rng.randrange(1800000))
                    for i in range(number))
                for left in positions:
                    for record in bam_record(chrom, left, index):
                        bam.write(record)
                    frag.write(f'chr{chrom:04}\t{left}\t{left + 150}\tBC{index % 8192:05}\t1\n'.encode())
                    index += 1
            bam.close()
            frag.close()
            for kind, path in paths.items():
                inputs[kind][role] = str(path)
        for kind in ('BAM', 'BAMPE', 'FRAG'):
            paths = inputs['BAM' if kind == 'BAMPE' else kind]
            cases.append({'name': f'{size}-{kind.lower()}', 'format': kind,
                          **paths, 'fragments_per_input': count,
                          'contigs': n_chrom, 'genome_size': n_chrom * 2000000,
                          'sha256': {role: sha256(path) for role, path in paths.items()}})
    manifest = {'schema': 1, 'generator': 'performance_fixtures.py/v1',
                'seed': 773, 'cases': cases}
    path = directory / 'manifest.json'
    path.write_text(json.dumps(manifest, indent=2) + '\n')
    return path
