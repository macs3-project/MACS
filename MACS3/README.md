# Directories

## Commands

The *Commands* directory  contains source codes for all subcommands.

1. callpeak_cmd.py
2. bdgpeakcall_cmd.py
3. bdgbroadcall_cmd.py
4. bdgcmp_cmd.py
5. bdgdiff_cmd.py
6. bdgopt_cmd.py
7. cmbreps_cmd.py
8. diffpeak_cmd.py
9. filterdup_cmd.py
10. pileup_cmd.py
11. predictd_cmd.py
12. randsample_cmd.py
13. refinepeak_cmd.py
14. callvar_cmd.py

## IO

The *IO* directory contains Python and Cython source codes for *read* and *write* data files.

1. BedGraphIO.pyx
2. OutputWriter.py
3. Parser.pyx
4. PeakIO.pyx

## Signal

The *Signal* directory contains Python, Cython and C source codes for data transformation and calculation.

1. BedGraph.pyx
2. CallPeakUnit.pyx
3. FixWidthTrack.pyx
4. PairedEndTrack.pyx
5. Pileup.pyx
6. Prob.pyx
7. ScoreTrack.pyx
8. Signal.pyx
9. cPosValCalculation.c, cPosValCalculation.h, cPosValCalculation.pxd
10. PeakModel.pyx
11. PeakDetect.pyx

## Utilities 

The *Utilities* directory contains source codes to validate command line options and to provide constants (including version number and buffer size).

1. Constants.py
2. OptValidator.py
