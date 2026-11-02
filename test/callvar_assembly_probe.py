"""Bounded native assembly probes, launched by test_callvar_assembly.py."""

from pathlib import Path
import sys

from MACS3.IO.PeakIO import PeakIO
from MACS3.Signal import RACollection as module
from MACS3.Signal.ReadAlignment import ReadAlignment


def reference(n=400):
    x = 20261003
    seq = bytearray()
    for _ in range(n):
        x = (1103515245 * x + 12345) % 2**31
        seq.append(b"ACGT"[(x >> 16) & 3])
    return bytes(seq)


REF = reference()


def read(start, seq, cigar=None, span=None):
    codes = [b"=ACMGRSVTWYHKDBN".index(c) for c in seq]
    if len(codes) % 2:
        codes.append(0)
    packed = bytes((codes[i] << 4) | codes[i + 1]
                   for i in range(0, len(codes), 2))
    return ReadAlignment(b"read", b"chr1", start,
                         start + (len(seq) if span is None else span), 0,
                         packed, bytes([40] * len(seq)),
                         cigar or (len(seq) << 4,), str(len(seq)).encode())


def collection(treatment, control=None):
    peak = PeakIO()
    peak.add(b"chr1", 0, 300)
    return module.RACollection(b"chr1", peak.get_data_from_chrom(b"chr1")[0],
                               treatment, control or [])


def main():
    assert Path(module.__file__).resolve() == Path(sys.argv[1]).resolve()
    assert Path(module.__file__).suffix != ".py"
    case = sys.argv[2]
    if case in ("short", "equal", "filtered", "ambiguous"):
        if case == "filtered":
            # No repeated k-mers: all reads are removed by fermi's filter.
            reads = [read(i, REF[i:i + 50]) for i in (0, 100, 200)]
        elif case == "ambiguous":
            reads = [read(i, b"N" * 50) for i in range(10)]
        else:
            reads = [read(i, REF[i:i + 25]) for i in range(10)]
        assert collection(reads).build_unitig_collection(
            25 if case == "equal" else 30) == 0
    elif case.startswith("invalid-"):
        overlap = int(case.split("-", 1)[1])
        c = collection([read(0, REF[:50])])
        try:
            c.build_unitig_collection(overlap)
        except ValueError as exc:
            assert "overlap" in str(exc).lower()
        else:
            raise AssertionError("invalid overlap accepted")
    elif case.startswith("mixed-") or case in ("minimum", "maximum"):
        long = [read(i, REF[i:i + 50]) for i in range(0, 249, 4)]
        short = [read(i, REF[i:i + 25]) for i in range(0, 249, 4)]
        if case == "mixed-treatment":
            c = collection(long + short)
        elif case == "mixed-control":
            c = collection(short, long)
        else:
            c = collection(long)
        overlap = 1 if case == "minimum" else 63 if case == "maximum" else 30
        uc = c.build_unitig_collection(overlap)
        if case == "maximum":
            assert uc == 0
        else:
            assert uc["URAs_list"]
            # Assembly can trim the ends when short reads are included.
            # The resulting mapped unitig must still match its reference.
            assert any(u["rpos"] - u["lpos"] >= 250
                       and u["seq"] == REF[u["lpos"]:u["rpos"]]
                       for u in uc["URAs_list"])
    elif case in ("bad-span", "splice"):
        r = (read(0, REF[:50], span=51) if case == "bad-span" else
             read(0, REF[:25] + REF[125:150],
                  cigar=(25 << 4, (100 << 4) | 3, 25 << 4), span=150))
        try:
            collection([r])
        except ValueError as exc:
            assert "reference span" in str(exc)
        else:
            raise AssertionError("inconsistent reference span accepted")
    else:
        raise AssertionError(case)


if __name__ == "__main__":
    main()
