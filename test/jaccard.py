#!/usr/bin/env python
# Time-stamp: <2024-02-12 15:21:42 Tao Liu>

import sys

from MACS3.Signal.Region import *


# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Calculate Jaccard Index of two peak files and return the index.\n\nneed 2 paras: %s <peak1> <peak2>\n" % sys.argv[0])
        sys.exit(1)

    #options.info("#  Read peak file 1...")
    peak1 = open( sys.argv[1] )
    regions1 = Regions()
    for l in peak1:
        fs = l.rstrip().split()
        regions1.add_loc( fs[0].encode(), int(fs[1]), int(fs[2]) )
    regions1.sort()

    #options.info("#  Read peak file 2...")
    peak2 = open( sys.argv[2] )
    regions2 = Regions()
    for l in peak2:
        fs = l.rstrip().split()
        regions2.add_loc( fs[0].encode(), int(fs[1]), int(fs[2]) )
    regions2.sort()

    tl1 = regions1.total_length()
    tl2 = regions2.total_length()
    inter = regions1.intersect(regions2)
    tl_inter = inter.total_length()

    jaccard_index = tl_inter / (tl1 + tl2 - tl_inter)
    print (jaccard_index)

if __name__ == '__main__':
    main()
