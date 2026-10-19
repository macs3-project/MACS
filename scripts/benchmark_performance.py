#!/usr/bin/env python3
"""Benchmark compiled MACS3 revisions without changing performance defaults.

Requires psutil in the runner only. Use prepare, run, then compare --help.
All target invocations run outside the source checkout with bounded timeouts.
"""

import argparse
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import sys
import time

from performance_fixtures import prepare, sha256

PROBE = Path(__file__).with_name('performance_probe.py').resolve()
SUFFIXES = ('_peaks.narrowPeak', '_peaks.broadPeak', '_peaks.gappedPeak',
            '_summits.bed', '_treat_pileup.bdg', '_control_lambda.bdg', '_peaks.xls')


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def artifacts(directory):
    """Hash stable bytes; only XLS comment metadata is excluded."""
    result = {}
    for path in sorted(Path(directory).iterdir()):
        if not path.name.endswith(SUFFIXES):
            continue
        data = path.read_bytes()
        if path.suffix == '.xls':
            data = b''.join(line for line in data.splitlines(keepends=True)
                            if not line.startswith(b'#'))
        result[path.name] = {'sha256': hashlib.sha256(data).hexdigest(),
                             'bytes': len(data)}
    required = {'bench_peaks.xls', 'bench_peaks.narrowPeak', 'bench_summits.bed',
                'bench_treat_pileup.bdg', 'bench_control_lambda.bdg'}
    if not required.issubset(result):
        raise RuntimeError(f'Missing scientific outputs: {required - result.keys()}')
    if not result['bench_peaks.narrowPeak']['bytes']:
        raise RuntimeError('Fixture produced no peaks')
    return result


def process_sample(pid, directory):
    """Sum tree RSS; deduplicate named/unlinked temp files by device+inode.

    Linux /proc exposes deleted open files. On macOS psutil.open_files cannot
    enumerate deleted files: report this gap, never call it a complete peak.
    Races with exit/unlink are expected. Permission errors are recorded.
    """
    import psutil
    files, errors = {}, set()
    rss = 0
    try:
        root = psutil.Process(pid)
        processes = [root, *root.children(recursive=True)]
    except psutil.NoSuchProcess:
        processes = []
    except psutil.AccessDenied:
        processes = []
        errors.add('process tree access denied')

    def record(path):
        try:
            st = path.stat()
            files[(st.st_dev, st.st_ino)] = st.st_size
        except FileNotFoundError:
            pass
        except PermissionError:
            errors.add('temp file access denied')

    for path in directory.rglob('*'):
        if path.is_file():
            record(path)
    prefix = str(directory) + os.sep
    for process in processes:
        try:
            rss += process.memory_info().rss
            if sys.platform == 'linux':
                fd_dir = Path(f'/proc/{process.pid}/fd')
                for fd in fd_dir.iterdir():
                    try:
                        if os.readlink(fd).startswith(prefix):
                            record(fd)
                    except FileNotFoundError:
                        pass
            else:
                for opened in process.open_files():
                    if opened.path.startswith(prefix):
                        record(Path(opened.path))
        except (psutil.NoSuchProcess, ProcessLookupError, FileNotFoundError):
            pass
        except (psutil.AccessDenied, PermissionError):
            errors.add('process metrics access denied')
    return rss, sum(files.values()), errors


def execute(command, cwd, env, log, timeout, interval, tempdir=None):
    start = time.perf_counter()
    peak_rss = peak_temp = samples = 0
    errors = set()
    timed_out = False
    with Path(log).open('w') as handle:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=handle,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while process.poll() is None:
                if tempdir is not None:
                    rss, used, issues = process_sample(process.pid, tempdir)
                    peak_rss, peak_temp = max(peak_rss, rss), max(peak_temp, used)
                    errors.update(issues)
                    samples += 1
                if time.perf_counter() - start > timeout:
                    timed_out = True
                    break
                time.sleep(interval)
        finally:
            # The session belongs exclusively to this invocation. Kill any
            # survivors even when its parent crashed before the timeout.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
    return {'status': process.returncode, 'timed_out': timed_out,
            'wall_seconds': time.perf_counter() - start,
            'sampled_tree_rss_bytes': peak_rss,
            'sampled_temp_bytes': peak_temp, 'samples': samples,
            'sampling_errors': sorted(errors), 'sample_interval_seconds': interval,
            'temp_includes_unlinked': sys.platform == 'linux',
            'command': command}


def target(args, env, directory, operation, spec, sample=False):
    filename = directory / f'{operation}.spec.json'
    result = directory / f'{operation}.result.json'
    dump(filename, {**spec, 'result': str(result), 'affinity': args.affinity,
                    'address_space_mib': args.address_space_mib})
    metrics = execute([args.python, str(PROBE), operation, str(filename)],
                      directory, env, directory / f'{operation}.log',
                      args.timeout, args.interval,
                      directory / 'tmp' if sample else None)
    dump(directory / f'{operation}.metrics.json', metrics)
    if metrics['status'] or not result.exists():
        raise RuntimeError(f'{operation} failed: {directory / (operation + ".log")}')
    return metrics, json.loads(result.read_text())


def run(args):
    if sys.platform not in ('darwin', 'linux'):
        raise ValueError('Supported measurement platforms: macOS and Linux')
    if (args.affinity or args.address_space_mib) and sys.platform != 'linux':
        raise ValueError('OS resource constraints require Linux; use a limited VM on macOS')
    manifest_path = Path(args.manifest).resolve()
    manifest = json.loads(manifest_path.read_text())
    cases = manifest['cases']
    if args.cases:
        selected = set(args.cases.split(','))
        cases = [case for case in cases if case['name'] in selected]
        if {case['name'] for case in cases} != selected:
            raise ValueError('Unknown case name')
    if not cases:
        raise ValueError('Empty workload matrix')
    for case in cases:
        for role in ('treat', 'control'):
            case[role] = str((manifest_path.parent / case[role]).resolve())
            if sha256(case[role]) != case['sha256'][role]:
                raise ValueError(f'Input checksum mismatch: {case[role]}')
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    if args.package_root:
        env['PYTHONPATH'] = str(Path(args.package_root).resolve())
    env['PYTHONHASHSEED'] = '0'
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        env[key] = '1'
    env.pop('MACS3_NUM_THREADS', None)
    metadata_dir = out / 'environment'
    metadata_dir.mkdir()
    _, metadata = target(args, env, metadata_dir, 'metadata', {})
    cli = str(Path(args.cli or metadata['cli']).resolve())
    metadata['cli'] = cli
    help_result = subprocess.run([args.python, cli, 'callpeak', '--help'],
                                 cwd=metadata_dir, env=env, capture_output=True,
                                 text=True, timeout=args.timeout)
    if help_result.returncode:
        raise RuntimeError(help_result.stderr)
    modern = '--threads' in help_result.stdout and '--process-backend' in help_result.stdout
    (metadata_dir / 'callpeak-help.txt').write_text(help_result.stdout)
    metadata.update({'label': args.label, 'revision': args.revision,
                     'profile': args.profile, 'supports_cpu_controls': modern,
                     'constraints': {'affinity': args.affinity,
                                     'address_space_mib': args.address_space_mib},
                     'measurement': 'tree RSS sums shared pages multiple times; sampled lower bounds; high-water RSS separately',
                     'manifest': manifest, 'runner_sha256': sha256(__file__),
                     'probe_sha256': sha256(PROBE)})
    dump(out / 'metadata.json', metadata)
    cache_dir = out / 'cache'
    cache_dir.mkdir()
    cache_metrics, cache_result = target(args, env, cache_dir, 'cache', {})
    dump(out / 'cache.json', {**cache_result, 'metrics': cache_metrics})
    policies = [(n, backend) for n in args.budgets for backend in ('serial', 'auto')]
    policies = [(n, b) for n, b in policies if n != 1 or b == 'serial'] if modern else [(None, None)]
    metadata['plan'] = [{'case': case['name'], 'budget': n, 'backend': b, 'repeat': repeat}
                        for repeat in range(1, args.repeats + 1)
                        for case in cases for n, b in policies]
    dump(out / 'metadata.json', metadata)
    runs = []
    for repeat in range(1, args.repeats + 1):
        # Reverse each repeat to distribute warm filesystem/order effects.
        order = [(case, n, b) for case in cases for n, b in policies]
        if repeat % 2 == 0:
            order.reverse()
        for case, budget, backend in order:
            name = f'{case["name"]}-{budget or "legacy"}-{backend or "default"}-{repeat}'
            directory = out / name
            directory.mkdir()
            tmp = directory / 'tmp'
            tmp.mkdir()
            output = directory / 'outputs'
            command_args = ['callpeak', '-t', case['treat'], '-c', case['control'],
                            '-f', case['format'], '-g', str(case['genome_size']),
                            '-n', 'bench', '--outdir', str(output), '--keep-dup', 'all',
                            '--nomodel', '--extsize', '150', '-B']
            if modern:
                command_args += ['--threads', str(budget), '--process-backend', backend]
            run_env = {**env, 'TMPDIR': str(tmp)}
            row = {'case': case['name'], 'budget': budget, 'backend': backend,
                   'repeat': repeat, 'profile': args.profile, 'label': args.label,
                   'revision': args.revision}
            print(f'{args.label}: {name}', flush=True)
            try:
                metrics, highwater = target(args, run_env, directory, 'cli',
                                            {'cli': cli, 'args': command_args}, sample=True)
                log = (directory / 'cli.log').read_text()
                row.update(metrics)
                row.update(highwater)
                row['execution_log'] = [line for line in log.splitlines()
                                        if re.search(r'Execution:|Peak process backend:|CPU limits:', line)]
                row['artifacts'] = artifacts(output)
                row['temp_remaining_bytes'] = sum(p.stat().st_size for p in tmp.rglob('*') if p.is_file())
                if row['temp_remaining_bytes']:
                    raise RuntimeError('temporary files retained after CLI exit')
                if row['sampling_errors']:
                    raise RuntimeError('resource sampling incomplete: ' + str(row['sampling_errors']))
            except (RuntimeError, subprocess.TimeoutExpired) as exc:
                row['error'] = str(exc)
            runs.append(row)
            dump(directory / 'run.json', row)
            dump(out / 'runs.json', runs)
    return 1 if any('error' in row for row in runs) else 0


def compare(args):
    groups, records, failures = defaultdict(list), [], []
    inputs = {}
    cache_limits = []
    for directory in map(Path, args.results):
        metadata = json.loads((directory / 'metadata.json').read_text())
        rows = json.loads((directory / 'runs.json').read_text())
        fields = ('case', 'budget', 'backend', 'repeat')
        actual = [tuple(row.get(field) for field in fields) for row in rows]
        planned = [tuple(row[field] for field in fields) for row in metadata['plan']]
        if not planned or sorted(actual, key=str) != sorted(planned, key=str):
            failures.append({'label': metadata['label'], 'error': 'missing, extra or duplicate runs'})
        for case in metadata['manifest']['cases']:
            identity = (case['format'], case['genome_size'], case['sha256'])
            if case['name'] in inputs and inputs[case['name']] != identity:
                failures.append({'case': case['name'], 'error': 'input definitions differ'})
            inputs[case['name']] = identity
        for row in rows:
            records.append(row)
            key = (row['label'], row['profile'], row['case'], row['budget'], row['backend'])
            groups[key].append(row)
            if 'error' in row:
                failures.append({'label': row['label'], 'case': row['case'], 'error': row['error']})
        cache = json.loads((directory / 'cache.json').read_text())
        phases = cache['phases']
        if (phases[-1]['live_table_bytes'] is None and
                metadata['revision'] == '170556de25475fc59b1e219e898c21b7eecd5dd8'):
            cache_limits.append({'label': metadata['label'], 'note': cache['instrumentation']})
        elif (phases[-1]['live_table_bytes'] != 64 * 1024**2 or
                phases[-1]['live_table_arrays'] != 2 or
                phases[-1]['allocations_to_date'] != phases[-2]['allocations_to_date'] or
                phases[-2]['allocations_to_date'] != phases[-3]['allocations_to_date']):
            failures.append({'label': metadata['label'], 'error': 'cache cap/reuse/lifetime gate failed'})
        if any(job['live_table_bytes'] != 64 * 1024**2 for job in cache.get('caller_jobs', [])):
            failures.append({'label': metadata['label'], 'error': 'repeated caller cache retention changed'})
    # Compare every repeat, budget, backend and revision within a profile/case.
    # OS resource profiles can also be compared explicitly by stable hashes.
    outputs = defaultdict(list)
    for row in records:
        if 'error' not in row:
            outputs[row['case']].append(row)
    equivalence = []
    for case, rows in sorted(outputs.items()):
        reference = rows[0]
        mismatches = [{'label': row['label'], 'budget': row['budget'],
                       'backend': row['backend'], 'repeat': row['repeat'],
                       'profile': row['profile']}
                      for row in rows if row['artifacts'] != reference['artifacts']]
        equivalence.append({'case': case, 'runs': len(rows),
                            'reference_label': reference['label'], 'mismatches': mismatches})
        if mismatches:
            failures.append({'case': case, 'error': 'scientific output bytes differ',
                             'mismatches': mismatches})
    lines = ['# Performance defaults benchmark', '',
             'RSS: sampled sum across live processes (shared pages counted repeatedly).',
             'Temp: sampled logical bytes; macOS excludes deleted open files. Both are lower bounds.', '',
             '| Revision | Profile | Case | Budget/backend | n | Wall median [min, max] s | Tree RSS MiB | Temp MiB |',
             '| --- | --- | --- | --- | ---: | --- | ---: | ---: |']
    for key, rows in sorted(groups.items(), key=lambda item: str(item[0])):
        good = [row for row in rows if 'error' not in row]
        if not good:
            continue
        label, profile, case, budget, backend = key
        wall = [row['wall_seconds'] for row in good]
        lines.append(f'| {label} | {profile} | {case} | {budget or "legacy"}/{backend or "default"} | '
                     f'{len(good)} | {statistics.median(wall):.3f} [{min(wall):.3f}, {max(wall):.3f}] | '
                     f'{statistics.median(r["sampled_tree_rss_bytes"] for r in good) / 1024**2:.2f} | '
                     f'{statistics.median(r["sampled_temp_bytes"] for r in good) / 1024**2:.2f} |')
    lines += ['', f'Scientific output/cache gates: {"FAIL" if failures else "PASS"}.', '',
              '```json', json.dumps({'equivalence': equivalence, 'failures': failures,
                                    'cache_observation_limits': cache_limits}, indent=2), '```', '']
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines))
    dump(path.with_suffix('.json'), {'equivalence': equivalence, 'failures': failures,
                                    'cache_observation_limits': cache_limits})
    return bool(failures)


def positive(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    fixture = commands.add_parser('prepare', help='Generate deterministic fixtures and SHA256 manifest')
    fixture.add_argument('--out', required=True)
    fixture.add_argument('--small', type=positive, default=6000)
    fixture.add_argument('--large', type=positive, default=500000)
    fixture.add_argument('--contigs', type=positive, default=256)
    runner = commands.add_parser('run', help='Measure one installed revision; never overwrite results')
    runner.add_argument('--manifest', required=True)
    runner.add_argument('--out', required=True)
    runner.add_argument('--label', required=True)
    runner.add_argument('--revision', required=True, help='Exact source SHA, recorded verbatim')
    runner.add_argument('--python', default=sys.executable)
    runner.add_argument('--cli', help='Selected revision bin/macs3; defaults to target interpreter scripts directory')
    runner.add_argument('--package-root', help='Explicit in-place package root; never inferred from checkout')
    runner.add_argument('--profile', default='host')
    runner.add_argument('--cases', help='Comma-separated manifest case names; default all')
    runner.add_argument('--budgets', type=lambda s: [positive(n) for n in s.split(',')], default=[1, 2, 4, 8])
    runner.add_argument('--repeats', type=positive, default=3)
    runner.add_argument('--timeout', type=positive, default=180)
    runner.add_argument('--interval', type=float, default=0.02)
    runner.add_argument('--affinity', type=lambda s: [int(n) for n in s.split(',')])
    runner.add_argument('--address-space-mib', type=positive, help='Linux per-process RLIMIT_AS, not a cgroup RSS limit')
    summary = commands.add_parser('compare', help='Fail on missing runs, output mismatch or cache lifetime regression')
    summary.add_argument('results', nargs='+')
    summary.add_argument('--out', required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        print(prepare(args.out, args.small, args.large, args.contigs))
        return 0
    if args.command == 'compare':
        return compare(args)
    if not 0 < args.interval <= 1:
        parser.error('--interval must be in (0, 1] seconds')
    return run(args)


if __name__ == '__main__':
    raise SystemExit(main())
