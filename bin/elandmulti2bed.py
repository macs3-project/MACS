#!/usr/bin/env python
# Time-stamp: <2009-11-08 20:03:20 Tao Liu>

import os
import sys
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Filter eland_multi.txt file and only keep the unique matches\nneed 2 paras: %s <eland_multi.txt> <bed output file>\n" % sys.argv[0])
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
    thisline = thisline.rstrip()
    if not thisline: return ("blank",None,None)

    if thisline[0] == "#": return ("comment line",None,None) # comment line is skipped
    thisfields = thisline.split()
    thistagname = thisfields[0]         # name of tag
    thistaglength = len(thisfields[1]) # length of tag

    if len(thisfields) < 4:
        return (None,None,None)
    else:
        thistaghits = sum(map(int,thisfields[2].split(':')))
        if thistaghits > 1:
            # multiple hits
            return (None,None,None)
        else:
            (name,pos) = thisfields[3].split(':')
            strand  = pos[-2]
            if strand == "F":
                fhd.write( "%s\t%d\t%d\t.\t.\t%s\n" % (name,
                                                       int(pos[:-2])-1,
                                                       int(pos[:-2])+thistaglength-1,
                                                       "+")
                           )
            elif strand == "R":
                fhd.write( "%s\t%d\t%d\t.\t.\t%s\n" % (name,
                                                       int(pos[:-2])-1,
                                                       int(pos[:-2])+thistaglength-1,
                                                       "-")
                           )
            else:
                raise StrandFormatError(thisline,strand)

class StrandFormatError(Exception):
    def __init__ (self, string, strand):
        self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)
        
    def __str__ (self):
        return repr(self.message)

if __name__ == '__main__':
    main()
