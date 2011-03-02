#!/usr/bin/env python
# Time-stamp: <2010-09-01 16:40:51 Tao Liu>

import os
import sys
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Convert ELAND_EXPORT file to BED\nneed 2 paras: %s <eland_export.txt> <bed output file>\n" % sys.argv[0])
        sys.exit(1)
    elandexport_fhd = file(sys.argv[1],"r")
    bed_fhd = file(sys.argv[2],"w")
    i = 0
    m = 0
    for thisline in elandexport_fhd:
        __fw_parse_line(thisline,bed_fhd)
        #i+=1
        #if i == 1000000:
        #    #m += 1
        #    #sys.stderr.write(" %d\n" % (m*1000000))
        #    #i=0

def __fw_parse_line (thisline, fhd ):
    if not thisline: return (None,None,None)
    if thisline.startswith("#") : return ("comment line",None,None) # comment line is skipped
    thisline = thisline.rstrip()
    if not thisline: return ("blank",None,None)
    
    thisfields = thisline.split("\t")

    if len(thisfields)>12 and thisfields[12]:
        thisname = ":".join(thisfields[0:6])
        thistaglength = len(thisfields[8])
        strand = thisfields[13]
        if strand == "F":
            fhd.write ( "%s\t%d\t%d\t%s\t.\t%s\n" % (thisfields[10],
                                                     int(thisfields[12])-1,
                                                     int(thisfields[12])+thistaglength-1,
                                                     thisname,
                                                     "+")
                        )
        elif strand == "R":
            fhd.write ( "%s\t%d\t%d\t%s\t.\t%s\n" % (thisfields[10],
                                                     int(thisfields[12])-1,
                                                     int(thisfields[12])+thistaglength-1,
                                                     thisname,
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
