# Time-stamp: <2012-03-05 15:11:06 Tao Liu>

"""Module for all MACS Parser classes for input.

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import logging
import struct
import gzip
from MACS14.Constants import *
from MACS14.IO.FeatIO import FWTrackII
# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "All Parser classes"

# ------------------------------------
# Misc functions
# ------------------------------------

def guess_parser ( fhd ):
    parser_dict = {"BED":BEDParser,
                   "ELAND":ELANDResultParser,
                   "ELANDMULTI":ELANDMultiParser,
                   "ELANDEXPORT":ELANDExportParser,
                   "SAM":SAMParser,
                   "BAM":BAMParser,
                   "BOWTIE":BowtieParser
                   }
    order_list = ("BAM",
                  "BED",
                  "ELAND",
                  "ELANDMULTI",
                  "ELANDEXPORT",
                  "SAM",
                  "BOWTIE",
                  )
    
    for f in order_list:
        p = parser_dict[f](fhd)
        s = p.sniff()
        if s:
            logging.info("Detected format is: %s" % (f) )
            return p
    raise Exception("Can't detect format!")

# ------------------------------------
# Classes
# ------------------------------------
class StrandFormatError(Exception):
    """Exception about strand format error.

    Example:
    raise StrandFormatError('Must be F or R','X')
    """
    def __init__ (self, string, strand):
        self.strand = strand
        self.string = string
        
    def __str__ (self):
        return repr( "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (self.string,self.strand) )
        
class GenericParser:
    """Generic Parser class.

    Inherit this to write your own parser.
    """
    def __init__ (self, fhd):
        self.fhd = fhd
        return

    def tsize(self):
        return

    def build_fwtrack (self):
        return 

    def __fw_parse_line (self, thisline ):
        return

    def sniff (self):
        try:
            t = self.tsize()
        except:
            self.fhd.seek(0)
            return False
        else:
            if t<=10 or t>=10000:
                self.fhd.seek(0)
                return False
            else:
                self.fhd.seek(0)
                return t

class BEDParser(GenericParser):
    """File Parser Class for tabular File.

    """
    def __init__ (self,fhd):
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split('\t')
            s += int(thisfields[2])-int(thisfields[1])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note: All locations will be merged (exclude the same
        location) then sorted after the track is built.

        If both_strand is True, it will store strand information in
        FWTrackII object.

        if do_merge is False, it will not merge the same location after
        the track is built.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        thisline = thisline.rstrip()
        if not thisline or thisline[:5]=="track" or thisline[:7]=="browser" or thisline[0]=="#": return ("comment line",None,None)

        thisfields = thisline.split('\t')
        chromname = thisfields[0]
        try:
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass

        if len(thisfields) < 6 : # default pos strand if no strand
                                 # info can be found
            return (chromname,
                    int(thisfields[1]),
                    0)
        else:
            if thisfields[5] == "+":
                return (chromname,
                        int(thisfields[1]),
                        0)
            elif thisfields[5] == "-":
                return (chromname,
                        int(thisfields[2]),
                        1)
            else:
                raise StrandFormatError(thisline,thisfields[5])

class ELANDResultParser(GenericParser):
    """File Parser Class for tabular File.

    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split('\t')
            s += len(thisfields[1])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self):
        """Build FWTrackII from all lines, return a FWTrackII object.

        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        #if thisline.startswith("#") or thisline.startswith("track") or thisline.startswith("browser"): return ("comment line",None,None) # comment line is skipped
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)

        thisfields = thisline.split('\t')
        thistaglength = len(thisfields[1])

        if len(thisfields) <= 6:
            return ("blank",None,None)

        try:
            chromname = thisfields[6]
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass

        if thisfields[2] == "U0" or thisfields[2]=="U1" or thisfields[2]=="U2":
            strand = thisfields[8]
            if strand == "F":
                return (chromname,
                        int(thisfields[7])-1,
                        0)
            elif strand == "R":
                return (chromname,
                        int(thisfields[7])+thistaglength-1,
                        1)
            else:
                raise StrandFormatError(thisline,strand)
        else:
            return (None,None,None)

class ELANDMultiParser(GenericParser):
    """File Parser Class for ELAND multi File.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields: 
    1. Sequence name 
    2. Sequence 
    3. Either NM, QC, RM (as described above) or the following: 
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches 
    found 
    5. Blank, if no matches found or if too many matches found, or the following: 
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction 
    starting at position 160322 with one error, one in the forward direction starting at 
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def tsize (self, strict = False):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split('\t')
            s += len(thisfields[1])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        if not thisline: return (None,None,None)
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)

        #if thisline[0] == "#": return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split('\t')
        thistagname = thisfields[0]        # name of tag
        thistaglength = len(thisfields[1]) # length of tag

        if len(thisfields) < 4:
            return (None,None,None)
        else:
            thistaghits = sum(map(int,thisfields[2].split(':')))
            if thistaghits > 1:
                # multiple hits
                return (None,None,None)
            else:
                (chromname,pos) = thisfields[3].split(':')

                try:
                    chromname = chromname[:chromname.rindex(".fa")]
                except ValueError:
                    pass
                
                strand  = pos[-2]
                if strand == "F":
                    return (chromname,
                            int(pos[:-2])-1,
                            0)
                elif strand == "R":
                    return (chromname,
                            int(pos[:-2])+thistaglength-1,
                            1)
                else:
                    raise StrandFormatError(thisline,strand)


class ELANDExportParser(GenericParser):
    """File Parser Class for ELAND Export File.

    """
    def __init__ (self,fhd):
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split("\t")
            s += len(thisfields[8])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        #if thisline.startswith("#") : return ("comment line",None,None) # comment line is skipped
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)
    
        thisfields = thisline.split("\t")

        if len(thisfields) > 12 and thisfields[12]:
            thisname = ":".join(thisfields[0:6])
            thistaglength = len(thisfields[8])
            strand = thisfields[13]
            if strand == "F":
                return (thisfields[10],int(thisfields[12])-1,0)
            elif strand == "R":
                return (thisfields[10],int(thisfields[12])+thistaglength-1,1)
            else:
                raise StrandFormatError(thisline,strand)
        else:
            return (None,None,None)

class PairEndELANDMultiParser(GenericParser):
    """File Parser Class for two ELAND multi Files for Pair-End sequencing.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields: 
    1. Sequence name 
    2. Sequence 
    3. Either NM, QC, RM (as described above) or the following: 
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches 
    found 
    5. Blank, if no matches found or if too many matches found, or the following: 
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction 
    starting at position 160322 with one error, one in the forward direction starting at 
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    def __init__ (self,lfhd,rfhd):
        self.lfhd = lfhd
        self.rfhd = rfhd        

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.lfhd.readline()
            thisline = thisline.rstrip()
            if not thisline: continue
            thisfields = thisline.split("\t")
            s += len(thisfields[1])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self, dist=200):
        """Build FWTrackII from all lines, return a FWTrackII object.

        lfhd: the filehandler for left tag file
        rfhd: the filehandler for right tag file
        dist: the best distance between two tags in a pair

        The score system for pairing two tags:

        score = abs(abs(rtag-ltag)-200)+error4lefttag+error4righttag

        the smaller score the better pairing. If the score for a
        pairing is bigger than 200, this pairing will be discarded.

        Note only the best pair is kept. If there are over two best
        pairings, this pair of left and right tags will be discarded.

        Note, the orders in left tag file and right tag file must
        match, i.e., the Nth left tag must has the same name as the
        Nth right tag.

        Note, remove comment lines beforehand.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        lnext = self.lfhd.next
        rnext = self.rfhd.next
        self.dist = dist
        try:
            while 1:
                lline = lnext()
                rline = rnext()
                (chromname,fpos,strand) = self.__fw_parse_line(lline,rline)

                i+=1
                if i == 1000000:
                    m += 1
                    logging.info(" %d" % (m*1000000))
                    i=0
                if not fpos or not chromname:
                    continue

                try:
                    chromname = chromname[:chromname.rindex(".fa")]
                except ValueError:
                    pass
                
                fwtrack.add_loc(chromname,fpos,strand)

        except StopIteration:
            pass
        return fwtrack
    
    def __fw_parse_line (self, leftline, rightline ):
        # >HWI-EAS275_5:4:100:340:1199/1	GTGCTGGTGGAGAGGGCAAACCACATTGACATGCT	2:1:0	chrI.fa:15061365F0,15068562F0,chrIV.fa:4783988R1
        # >HWI-EAS275_5:4:100:340:1199/2	GGTGGTGTGTCCCCCTCTCCACCAGCACTGCGGCT	3:0:0	chrI.fa:15061451R0,15068648R0,15071742R0

        leftfields = leftline.split('\t')
        lefttaglength = len(leftfields[1]) # length of tag
        rightfields = rightline.split('\t')
        righttaglength = len(rightfields[1]) # length of tag

        if len(rightfields) < 4 or len(leftfields) < 4:
            # one of the tag cann't be mapped to genome
            return (None,None,None)
        else:
            lefthits = self.__parse_line_to_dict(leftfields[3])
            righthits = self.__parse_line_to_dict(rightfields[3])            
            parings = []

            for seqname in lefthits.keys():
                if not righthits.has_key(seqname):
                    continue
                else:
                    leftpses = lefthits[seqname] # pse=position+strand+error
                    rightpses = righthits[seqname]
                    for (lp,ls,le) in leftpses:
                        for (rp,rs,re) in rightpses:
                            # try to pair them
                            if ls == 'F':
                                if rs == 'R':
                                    score = abs(abs(rp-lp)-self.dist)+le+re
                                    if score < 200:
                                        #parings.append((score,seqname,int((lp+rp)/2),0) )
                                        parings.append((score,seqname,lp,0))
                                else:
                                    # strands don't match
                                    continue
                            else:
                                if rs == 'F':
                                    score = abs(abs(rp-lp)-self.dist)+le+re
                                    if score < 200:
                                        #parings.append((score,seqname,int((lp+rp)/2),1) )
                                        parings.append((score,seqname,lp,1))
                                else:
                                    # strands don't match
                                    continue
            if not parings:
                return (None,None,None)
            parings.sort()
            if len(parings)>1 and parings[0][0] == parings[1][0]:
                # >2 best paring, reject!
                return (None,None,None)
            else:
                return parings[0][1:]                                
                    
    def __parse_line_to_dict ( self, linestr ):
        items = linestr.split(',')
        hits = {}
        for item in items:
            if item.find(':') != -1:
                # a seqname section
                (n,pse) = item.split(":") # pse=position+strand+error
                try:
                    n = n[:n.rindex(".fa")]
                except ValueError:
                    pass

                hits[n]=[]
                try:
                    sindex = pse.rindex('F')
                except ValueError:
                    sindex = pse.rindex('R')
                p = int(pse[:sindex])
                s = pse[sindex]
                e = int(pse[sindex+1:])
                hits[n].append((p,s,e))
            else:
                # only pse section
                try:
                    sindex = pse.rindex('F')
                except ValueError:
                    sindex = pse.rindex('R')
                p = int(pse[:sindex])
                s = pse[sindex]
                e = int(pse[sindex+1:])
                hits[n].append((p,s,e))
        return hits

### Contributed by Davide, modified by Tao
class SAMParser(GenericParser):
    """File Parser Class for SAM File.

    Each line of the output file contains at least: 
    1. Sequence name 
    2. Bitwise flag
    3. Reference name
    4. 1-based leftmost position fo clipped alignment
    5. Mapping quality
    6. CIGAR string
    7. Mate Reference Name
    8. 1-based leftmost Mate Position
    9. Inferred insert size
    10. Query sequence on the same strand as the reference
    11. Query quality
    
    The bitwise flag is made like this:
    dec	meaning
    ---	-------
    1	paired read
    2	proper pair
    4	query unmapped
    8	mate unmapped
    16	strand of the query (1 -> reverse)
    32	strand of the mate
    64	first read in pair
    128	second read in pair
    256	alignment is not primary
    512	does not pass quality check
    1024	PCR or optical duplicate
    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split("\t")
            s += len(thisfields[9])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)
        if thisline[0]=="@": return ("comment line",None,None) # header line started with '@' is skipped
        thisfields = thisline.split('\t')
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
            # From Benjamin Schiller https://github.com/benjschiller
            if bwflag & 128:
                # this is not the first read in a pair
                return (None, -1, None)
            # end of the patch
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            thisstrand = 1
            thisstart = int(thisfields[3]) - 1 + len(thisfields[9])	#reverse strand should be shifted len(query) bp 
        else:
            thisstrand = 0
            thisstart = int(thisfields[3]) - 1	

        try:
            thisref = thisref[:thisref.rindex(".fa")]
        except ValueError:
            pass
        return (thisref, thisstart, thisstrand)

class BAMParser(GenericParser):
    """File Parser Class for BAM File.

    File is gzip-compatible and binary.
    Information available is the same that is in SAM format.
    
    The bitwise flag is made like this:
    dec	meaning
    ---	-------
    1	paired read
    2	proper pair
    4	query unmapped
    8	mate unmapped
    16	strand of the query (1 -> reverse)
    32	strand of the mate
    64	first read in pair
    128	second read in pair
    256	alignment is not primary
    512	does not pass quality check
    1024	PCR or optical duplicate
    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def sniff(self):
        if self.fhd.read(3) == "BAM":
            return True
        else:
            return False
            
    def tsize(self):
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek(4)
        header_len =  struct.unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = struct.unpack('<i', fread(4))[0]
        #print nc
        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fread(4))[0]
            # jump over chromosome size, we don't need it
            fread(nlength)
            fseek(ftell() + 4)
        s = 0
        n = 0
        while n<10:
            #print n
            entrylength = struct.unpack('<i', fread(4))[0]
            #print entrylength
            data = fread(entrylength)
            #print n
            #(chrid,fpos,strand) = self.__fw_binary_parse(fread(entrylength))
            s += struct.unpack('<i', data[16:20])[0]
            #print n
            n += 1
            #print s,n
        fseek(0)
        logging.info("tag size: %d" % int(s/n))
        return int(s/n)

    def build_fwtrack (self):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        references = []
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek(4)
        header_len =  struct.unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = struct.unpack('<i', fread(4))[0]
        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fread(4))[0]
            references.append(fread(nlength)[:-1])
            # jump over chromosome size, we don't need it
            fseek(ftell() + 4)
        
        while 1:
            try:
                entrylength = struct.unpack('<i', fread(4))[0]
            except struct.error:
                break
            (chrid,fpos,strand) = self.__fw_binary_parse(fread(entrylength))                    
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if fpos >= 0:
                fwtrack.add_loc(references[chrid],fpos,strand)
        self.fhd.close()
        return fwtrack
    
    def __fw_binary_parse (self, data ):
        # we skip lot of the available information in data (i.e. tag name, quality etc etc)
        if not data: return (None,-1,None)

        thisref = struct.unpack('<i', data[0:4])[0]
        thisstart = struct.unpack('<i', data[4:8])[0]
        (cigar, bwflag) = struct.unpack('<hh', data[12:16])
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return (None, -1, None)       #unmapped sequence or bad sequence
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return (None, -1, None)   # not a proper pair
            if bwflag & 8:
                return (None, -1, None)   # the mate is unmapped
            # From Benjamin Schiller https://github.com/benjschiller
            if bwflag & 128:
                # this is not the first read in a pair
                return (None, -1, None)
            # end of the patch
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            thisstrand = 1
            thisstart = thisstart + struct.unpack('<i', data[16:20])[0]	#reverse strand should be shifted len(query) bp 
        else:
            thisstrand = 0

        return (thisref, thisstart, thisstrand)

### End ###

class BowtieParser(GenericParser):
    """File Parser Class for map files from Bowtie or MAQ's maqview
    program.

    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split('\t')
            s += len(thisfields[4])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note: All locations will be merged (exclude the same
        location) then sorted after the track is built.

        If both_strand is True, it will store strand information in
        FWTrackII object.

        if do_merge is False, it will not merge the same location after
        the track is built.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        """
        The following definition comes from bowtie website:
        
        The bowtie aligner outputs each alignment on a separate
        line. Each line is a collection of 8 fields separated by tabs;
        from left to right, the fields are:

        1. Name of read that aligned

        2. Orientation of read in the alignment, - for reverse
        complement, + otherwise

        3. Name of reference sequence where alignment occurs, or
        ordinal ID if no name was provided

        4. 0-based offset into the forward reference strand where
        leftmost character of the alignment occurs

        5. Read sequence (reverse-complemented if orientation is -)

        6. ASCII-encoded read qualities (reversed if orientation is
        -). The encoded quality values are on the Phred scale and the
        encoding is ASCII-offset by 33 (ASCII char !).

        7. Number of other instances where the same read aligns
        against the same reference characters as were aligned against
        in this alignment. This is not the number of other places the
        read aligns with the same number of mismatches. The number in
        this column is generally not a good proxy for that number
        (e.g., the number in this column may be '0' while the number
        of other alignments with the same number of mismatches might
        be large). This column was previously described as"Reserved".

        8. Comma-separated list of mismatch descriptors. If there are
        no mismatches in the alignment, this field is empty. A single
        descriptor has the format offset:reference-base>read-base. The
        offset is expressed as a 0-based offset from the high-quality
        (5') end of the read.

        """
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)
        if thisline[0]=="#": return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split('\t')                      # I hope it will never bring me more trouble

        chromname = thisfields[2]
        try:
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass

            if thisfields[1] == "+":
                return (chromname,
                        int(thisfields[3]),
                        0)
            elif thisfields[1] == "-":
                return (chromname,
                        int(thisfields[3])+len(thisfields[4]),
                        1)
            else:
                raise StrandFormatError(thisline,thisfields[1])

