#!/usr/bin/env python
# Time-stamp: <2010-06-07 17:34:30 Tao Liu>

import os
import sys
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Convert SAM file to BED\nneed 2 paras: %s <sam> <bed output file>\n" % sys.argv[0])
        sys.exit(1)
    sam_fhd = file(sys.argv[1],"r")
    bed_fhd = file(sys.argv[2],"w")
    i = 0
    m = 0
    for thisline in sam_fhd:
        __fw_parse_line(thisline,bed_fhd)

def __fw_parse_line ( thisline,bed_fhd ):
    thisline = thisline.rstrip()
    if not thisline: return ("blank",None,None)
    if thisline[0]=="@": return ("comment line",None,None) # header line started with '@' is skipped
    thisfields = thisline.split()
    thistagname = thisfields[0]         # name of tag
    thisref = thisfields[2]
    bwflag = int(thisfields[1])
    if bwflag & 4 or bwflag & 512 or bwflag & 1024:
        return (None, None, None)       #unmapped sequence or bad sequence
    if bwflag & 1:
        # paired read. We should only keep sequence if the mate is mapped
        # and if this is the left mate, all is within  the flag! 
        if not bwflag & 2:
            return (None, None, None)   # not a proper pair
        if bwflag & 8:
            return (None, None, None)   # the mate is unmapped
        p1pos = int(thisfields[3]) - 1
        p2pos = int(thisfields[7]) - 1
        if p1pos > p2pos:
            # this pair is the farthest one, skip it
            return (None, None, None)
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
    try:
        thisref = thisref[:thisref.rindex(".fa")]
    except ValueError:
        pass

    if bwflag & 16:
        bed_fhd.write( "%s\t%d\t%d\t%s\t%d\t%s\n" % (thisref,int(thisfields[3])-1,int(thisfields[3])-1+len(thisfields[9]),thisfields[0],0,"-") )
    else:
        bed_fhd.write( "%s\t%d\t%d\t%s\t%d\t%s\n" % (thisref,int(thisfields[3])-1,int(thisfields[3])-1+len(thisfields[9]),thisfields[0],0,"+") )


class StrandFormatError(Exception):
    def __init__ (self, string, strand):
        self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

    def __str__ (self):
        return repr(self.message)

if __name__ == '__main__':
    main()
