# CPU budgets and process backends

MACS3 can limit the CPU resources used by supported commands. Use
`--threads N` to set a budget for one invocation:

```sh
macs3 callpeak -t treatment.bam -c control.bam -g hs --threads 4
macs3 callpeak -t treatment.bam -g hs --threads 1
```

The budget includes the main process or thread. A value of `1` runs MACS3
serially. The option is available for `callpeak`, `predictd`, `pileup`,
`filterdup`, `randsample`, `refinepeak`, `hmmratac`, and `callvar`.

When `--threads` is not provided, MACS3 uses `MACS3_NUM_THREADS` if set to a
positive integer. Otherwise, it requests eight CPUs, limited to the resources
available to the process. `callvar` keeps its historical default of one CPU.
An explicit `--threads` value takes precedence over the environment variable.

`callpeak` and `callvar` also accept `--process-backend {auto,serial}`. The
default, `auto`, uses process workers where supported; `serial` disables
process workers while allowing other work to use the requested thread budget.
Use `--threads 1` for fully serial execution.

## Python API

Python callers can pass an execution policy to parsers, tracks, peak models,
and peak callers. Resolve a policy once and share it across objects belonging
to the same invocation:

```python
from MACS3.Utilities.Execution import resolve_execution_policy
from MACS3.IO.Parser import BAMParser

policy = resolve_execution_policy(threads=4, process_backend="serial")
parser = BAMParser("treatment.bam", execution_policy=policy)
track = parser.build_fwtrack()
track.finalize()
```

The policy can also be supplied as `threads=` where supported. See the
[Python API reference](api/index) for public classes and functions.
