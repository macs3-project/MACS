#!/usr/bin/env bash

## pyprofile: profile Python script

## $Revision$
## Copyright 2011 Michael M. Hoffman <mmh1@uw.edu>
## Modified by Tao Liu: Correctly get the command name.

set -o nounset
set -o pipefail
set -o errexit

if [ $# -lt 2 ]; then
    echo usage: "$0" PROF CMDLINE...
    exit 2
fi

PROF=$1
MAIN=`basename $2`
CMD=( $@ )
CMD=${CMD[@]: 1}

echo $PROF
echo $MAIN
echo $CMD

python -m cProfile -o $PROF $CMD
./pyprofile_stat.py $PROF tottime > ${PROF/.prof/}.tottime
./pyprofile_stat.py $PROF calls > ${PROF/.prof/}.calls
./pyprofile_stat.py $PROF cumulative > ${PROF/.prof/}.cumulative
