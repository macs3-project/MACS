#!/usr/bin/env python

from MACS2.IO.Parser import BEDParser

#fhd = gzip.open("peakcalling/ChIP_0.1.bam","r")

a_parser = BEDParser("ChIP_0.1.bed.gz")

a = parser.build_fwtrack()

b_parser = BEDParser("Control_0.1.bed.gz")

b = parser.build_fwtrack()

