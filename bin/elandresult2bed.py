#!/usr/bin/env python
# Time-stamp: <2009-11-08 20:03:04 Tao Liu>

import os
import sys
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Filter eland_result.txt file and only keep the unique matches\nneed 2 paras: %s <eland_result.txt> <bed output file>\n" % sys.argv[0])
        sys.exit(1)
    eland_multi_fhd = file(sys.argv[1],"r")
    bed_fhd = file(sys.argv[2],"w")
    i = 0
    m = 0
    for thisline in eland_multi_fhd:
        __fw_parse_line(thisline,bed_fhd)
        i+=1
        if i == 1000000:
            m += 1
            sys.stderr.write(" %d\n" % (m*1000000))
            i=0

def __fw_parse_line (thisline, fhd ):
    if not thisline: return (None,None,None)
    if thisline.startswith("#") or thisline.startswith("track") or thisline.startswith("browser"): return ("comment line",None,None) # comment line is skipped
    thisline = thisline.rstrip()
    if not thisline: return ("blank",None,None)
    
    thisfields = thisline.split()
    thistaglength = len(thisfields[1])
    
    if thisfields[2] == "U0" or thisfields[2]=="U1" or thisfields[2]=="U2":
        strand = thisfields[8]
        if strand == "F":
            fhd.write ( "%s\t%d\t%d\t.\t.\t%s\n" % (thisfields[6],
                                                    int(thisfields[7])-1,
                                                    int(thisfields[7])+thistaglength-1,
                                                    "+")
                        )
        elif strand == "R":
            fhd.write ( "%s\t%d\t%d\t.\t.\t%s\n" % (thisfields[6],
                                                    int(thisfields[7])-1,
                                                    int(thisfields[7])+thistaglength-1,
                                                    "-")
                        )
        else:
            raise StrandFormatError(thisline,strand)
    else:
        return (None,None,None)

class StrandFormatError(Exception):
    def __init__ (self, string, strand):
        self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)
        
    def __str__ (self):
        return repr(self.message)

if __name__ == '__main__':
    main()
