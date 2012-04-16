#!/usr/bin/env python

from MACS2.IO.cParser import BEDParser, BAMParser
import gzip

#fhd = gzip.open("peakcalling/ChIP_0.1.bam","r")

parser = BAMParser("peakcalling/ChIP_0.1.bam")

parser.build_fwtrack()

