#!/usr/bin/env python
# Time-stamp: <2019-10-30 16:07:15 taoliu>

import sys
import pstats
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("need 2 paras: %s <prof file> <time|calls|cumulative|...>\n" % sys.argv[0])
        sys.exit(1)

    p = pstats.Stats(sys.argv[1])
    sort_term = sys.argv[2]

    p.strip_dirs().sort_stats(sort_term).print_stats(20)


if __name__ == '__main__':
    main()
    
