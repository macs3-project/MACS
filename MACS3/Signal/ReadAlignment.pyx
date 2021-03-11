# cython: language_level=3
# cython: profile=True
# Time-stamp: <2021-03-10 16:21:51 Tao Liu>

"""Module for SAPPER ReadAlignment class

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from cpython cimport bool

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)
    void *malloc(size_t size)
    void *calloc(size_t n, size_t size)
    void free(void *ptr)
    int strcmp(char *a, char *b)
    char * strcpy(char *a, char *b)
    long atol(char *bytes)
    int atoi(char *bytes)

# ------------------------------------
# constants
# ------------------------------------

__BAMDNACODE__ = b"=ACMGRSVTWYHKDBN"
__CIGARCODE__ = "MIDNSHP=X"
__DNACOMPLEMENT__ = b'\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TBGDEFCHIJKLMNOPQRSAUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff' # A trans table to convert A to T, C to G, G to C, and T to A.

# -- CIGAR CODE --
#OP BAM  Description
#M  0    alignment match (can be a sequence match or mismatch) insertion to the reference
#I  1    insertion to the reference
#D  2    deletion from the reference
#N  3    skipped region from the reference
#S  4    soft clipping (clipped sequences present in SEQ)
#H  5    hard clipping (clipped sequences NOT present in SEQ)
#P  6    padding (silent deletion from padded reference)
#=  7    sequence match
#X  8    sequence mismatch
# -- -- -- -- -- --

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
cdef class ReadAlignment:
    cdef:
        bytes readname
        bytes chrom
        int lpos
        int rpos
        int strand              # strand information. 0 means forward strand, 1 means reverse strand.
        bytes binaryseq
        bytes binaryqual
        int l                   # length of read
        tuple cigar             # each item contains op_l|op
        bytes MD
        int n_edits             # number of edits; higher the number,
                                # more differences with reference
        bytes SEQ               # sequence of read regarding to + strand
        bytes QUAL              # quality of read regarding to + strand

    def __init__ ( self,
                   bytes readname,
                   bytes chrom, int lpos, int rpos,
                   int strand,
                   bytes binaryseq, 
                   bytes binaryqual,
                   tuple cigar,
                   bytes MD ):
        self.readname = readname
        self.chrom = chrom
        self.lpos = lpos
        self.rpos = rpos
        self.strand = strand
        self.binaryseq = binaryseq
        self.binaryqual = binaryqual
        self.l = len( binaryqual )
        self.cigar = cigar
        self.MD = MD
        self.n_edits = self.get_n_edits()
        (self.SEQ, self.QUAL) = self.__get_SEQ_QUAL()

    cdef int get_n_edits( self ):
        """The number is from self.cigar and self.MD.

        """
        cdef:
            int n_edits
            int i, cigar_op, cigar_op_l
            char c
        
        n_edits = 0
        for i in self.cigar:    # only count insertion or softclip
            cigar_op = i & 15
            cigar_op_l = i >> 4
            if cigar_op in [ 1, 4 ]:    # count Insertion or Softclip
                n_edits += cigar_op_l
        
        for c in self.MD:
            if (c > 64 and c < 91): # either deletion in query or mismatch
                n_edits += 1
        return n_edits

    def __str__ ( self ):
        c = self.chrom.decode()
        n = self.readname.decode()
        if self.strand:
            s = "-"
        else:
            s = "+"
        return f"{c}\t{self.lpos}\t{self.rpos}\t{n}\t{self.l}\t{s}"
    
    def __getitem__ ( self, keyname ):
        if keyname == "readname":
            return self.readname
        elif keyname == "chrom":
            return self.chrom
        elif keyname == "lpos":
            return self.lpos
        elif keyname == "rpos":
            return self.rpos
        elif keyname == "strand":
            return self.strand
        elif keyname == "SEQ":
            return self.SEQ
        elif keyname == "QUAL":
            return self.QUAL
        elif keyname == "n_edits":
            return self.n_edits
        elif keyname == "binaryseq":
            return self.binaryseq
        elif keyname == "binaryqual":
            return self.binaryqual
        elif keyname == "l":
            return self.l
        elif keyname == "cigar":
            return self.cigar
        elif keyname == "MD":
            return self.MD
        else:
            raise KeyError("No such key", keyname)

    def __getstate__ ( self ):
        return ( self.readname, self.chrom, self.lpos, self.rpos, self.strand, self.binaryseq, self.binaryqual, self.l, self.cigar, self.MD, self.n_edits, self.SEQ, self.QUAL )

    def __setstate__ ( self, state ):
        ( self.readname, self.chrom, self.lpos, self.rpos, self.strand, self.binaryseq, self.binaryqual, self.l, self.cigar, self.MD, self.n_edits, self.SEQ, self.QUAL ) = state

    # cpdef bytearray get_SEQ ( self ):
    #     """Convert binary seq to ascii seq.

    #     Rule: for each byte, 1st base in the highest 4bit; 2nd in the lowest 4bit. "=ACMGRSVTWYHKDBN" -> [0,15]

    #     Note: In BAM, if a sequence is mapped to reverse strand, the
    #     reverse complement seq is written in SEQ field. So the return
    #     value of this function will not be the original one if the
    #     read is mapped to - strand.
    #     """
    #     cdef:
    #         char c
    #         bytearray seq

    #     seq = bytearray(b"")
    #     for c in self.binaryseq:
    #         # high
    #         seq.append( __BAMDNACODE__[c >> 4 & 15] )
    #         # low
    #         seq.append( __BAMDNACODE__[c & 15] )
    #     if seq[-1] == b"=":
    #         # trim the last '=' if it exists
    #         seq = seq[:-1]
    #     return seq

    cdef tuple __get_SEQ_QUAL ( self ):
        """Get tuple of (SEQ, QUAL).

        Rule: for each byte, 1st base in the highest 4bit; 2nd in the lowest 4bit. "=ACMGRSVTWYHKDBN" -> [0,15]

        Note: In BAM, if a sequence is mapped to reverse strand, the
        reverse complement seq is written in SEQ field. So the return
        value of this function will not be the original one if the
        read is mapped to - strand. If you need to original one, do
        reversecomp for SEQ and reverse QUAL.
        """
        cdef:
            int i
            char c
            bytearray seq
            bytearray qual

        seq = bytearray(b"")
        qual = bytearray(b"")

        for i in range( len(self.binaryseq) ):
            c = self.binaryseq[ i ]
            # high
            seq.append( __BAMDNACODE__[c >> 4 & 15] )
            # low
            seq.append( __BAMDNACODE__[c & 15] )

        for i in range( len( self.binaryqual ) ):
            # qual is the -10log10 p or phred score. 
            qual.append( self.binaryqual[i] )
            
        if seq[-1] == b"=":
            # trim the last '=' if it exists
            seq = seq[:-1]
        assert len( seq ) == len( qual ), Exception("Lengths of seq and qual are not consistent!")

        # Example on how to get original SEQ and QUAL:
        #if self.strand:
        #    seq.reverse()
        #    #compliment
        #    seq = seq.translate( __DNACOMPLEMENT__ )
        #    qual.reverse()

        return ( bytes(seq), bytes(qual) )

    
    cpdef bytes get_FASTQ ( self ):
        """Get FASTQ format text.

        """
        cdef:
            bytes seq
            bytearray qual

        seq = self.SEQ
        qual = bytearray(self.QUAL)

        for i in range( len( self.QUAL ) ):
            # qual is the -10log10 p or phred score, to make FASTQ, we have to add 33
            qual[ i ] += 33
        
        # reverse while necessary
        if self.strand:
            seq = self.SEQ[::-1]
            #compliment
            seq = seq.translate( __DNACOMPLEMENT__ )
            qual = qual[::-1]
        else:
            seq = self.SEQ

        return b"@" + self.readname + b"\n" + seq + b"\n+\n" + qual + b"\n"

    cpdef bytearray get_REFSEQ ( self ):
        """Fetch reference sequence, using self.MD and self.cigar
        """
        cdef:
            char c
            bytearray seq, refseq
            int i, cigar_op, cigar_op_l
            bytearray MD_op
            int ind
            bool flag_del       # flag for deletion event in query

        seq = bytearray(self.SEQ)    # we start with read seq then make modifications

        # 2-step proces
        # First step: use CIGAR to edit SEQ to remove S (softclip) and I (insert)
        # __CIGARCODE__ = "MIDNSHP=X"
        # let ind be the index in SEQ
        ind = 0
        for i in self.cigar:
            cigar_op = i & 15
            cigar_op_l = i >> 4
            if cigar_op in [2, 5, 6]:     # do nothing for Deletion (we will
                                          # put sequence back in step 2),
                                          # Hardclip and Padding
                pass
            elif cigar_op in [0, 7, 8]:   # M = X alignment match (match or
                                          # mismatch)
                # do nothing and move ind
                ind += cigar_op_l
            elif cigar_op in [ 1, 4 ]:    # Remove for Insertion or Softclip
                seq[ ind : ind + cigar_op_l ] = b''

        # now the seq should be at the same length as rpos-lpos 

        # Second step: use MD string to edit SEQ to put back 'deleted
        # seqs' and modify mismatches

        # let ind be the index in SEQ again, from 0
        ind = 0
        MD_op = bytearray(b'')
        flag_del = False
        for c in self.MD:
            if c < 58 and c > 47:
                # means Match
                flag_del = False
                MD_op.append(c)
            elif (c > 64 and c < 91) and not flag_del:
                # An alphabet means Mismatch, Note, if MD is made
                # right, a mismatch should only be 1 letter surrounded
                # by digits.
                ind += int(MD_op)
                seq[ ind ] = c
                ind += 1
                # reset MD_op
                MD_op = bytearray(b'')
            elif (c > 64 and c < 91) and flag_del:
                seq[ ind:ind ] = [c,]
                ind += 1
            elif c == 94:
                # means Deletion in query. Now, insert a sequnce into
                # SEQ
                flag_del = True
                ind += int(MD_op)
                # reset MD_op
                MD_op = bytearray(b'')
            else:
                raise Exception("Don't understand this operator in MD: %c" % c)
            #print( seq.decode() )

        return seq
    
    cpdef get_base_by_ref_pos ( self, long ref_pos ):
        """Get base by ref position.

        """
        cdef:
           int relative_pos, p
        assert self.lpos <= ref_pos and self.rpos > ref_pos, Exception("Given position out of alignment location")
        relative_pos = ref_pos - self.lpos
        p = self.relative_ref_pos_to_relative_query_pos( relative_pos ) 

        if p == -1:             # located in a region deleted in query
            return None
        else:
            return __BAMDNACODE__[ (self.binaryseq[p//2] >> ((1-p%2)*4) ) & 15 ]

    cpdef get_bq_by_ref_pos ( self, long ref_pos ):
        """Get base quality by ref position. Base quality is in Phred scale.

        Returned value is the raw Phred-scaled base quality.

        """
        cdef:
           int relative_pos, p
        assert self.lpos <= ref_pos and self.rpos > ref_pos, Exception("Given position out of alignment location")
        relative_pos = ref_pos - self.lpos
        p = self.relative_ref_pos_to_relative_query_pos( relative_pos ) 

        if p == -1:             # located in a region deleted in query
            return None
        else:
            return self.binaryqual[p]

    cpdef tuple get_base_bq_by_ref_pos ( self, long ref_pos ):
        """Get base and base quality by ref position. Base quality is in Phred scale.

        Returned bq is the raw Phred-scaled base quality.
        """
        cdef:
           int relative_pos, p
        assert self.lpos <= ref_pos and self.rpos > ref_pos, Exception("Given position out of alignment location")
        relative_pos = ref_pos - self.lpos
        p = self.relative_ref_pos_to_relative_query_pos( relative_pos ) 

        if p == -1:             # located in a region deleted in query
            return None
        else:
            return ( __BAMDNACODE__[ (self.binaryseq[p//2] >> ((1-p%2)*4) ) & 15 ], self.binaryqual[p] )

    cpdef tuple get_variant_bq_by_ref_pos ( self, long ref_pos ):
        """Get any variants (different with reference) and base quality by ref position. 

        variants will be 

        1) =, if identical

        2) A/T/C/G, if SNV

        3) -, if the reference base is deleted, in this case, bq will
        be the highest possible bq, which is 93.

        4) ^<A/T/C/G>+, if there is an insertion at the location

        Base quality is the raw Phred-scaled base quality.

        """
        cdef:
           int i, m, n
           int res, p, op, op_l
           int pos
           bool tip
           bytearray refseq
           bytes p_refseq, p_seq
           bytearray seq_array
           bytearray bq_array

        assert self.lpos <= ref_pos and self.rpos > ref_pos, Exception("Given position out of alignment location")

        res = ref_pos - self.lpos          # residue
        p = 0
        
        refseq = self.get_REFSEQ()
        p_refseq =  refseq[ res ]
        # -- CIGAR CODE --
        #OP BAM  Description
        #M  0    alignment match (can be a sequence match or mismatch) insertion to the reference
        #I  1    insertion to the reference
        #D  2    deletion from the reference
        #N  3    skipped region from the reference
        #S  4    soft clipping (clipped sequences present in SEQ)
        #H  5    hard clipping (clipped sequences NOT present in SEQ)
        #P  6    padding (silent deletion from padded reference)
        #=  7    sequence match
        #X  8    sequence mismatch

        seq_array = bytearray( b'' )
        bq_array = bytearray( b'' )

        for m in range( len(self.cigar) ):
            i = self.cigar[ m ]
            op = i & 15
            op_l = i >> 4
            if op in [0, 7, 8]:         # M = X alignment match (match or mismatch)
                if res < op_l - 1:
                    # in the range of a CIGAR operator
                    p += res
                    # find the position, now get the ref
                    pos = p
                    seq_array.append( __BAMDNACODE__[ (self.binaryseq[ p//2 ] >> ( (1-p%2)*4 ) ) & 15 ] )
                    bq_array.append( self.binaryqual[ p ] )
                    break
                elif res == op_l - 1:
                    p += res
                    pos = p
                    seq_array.append( __BAMDNACODE__[ (self.binaryseq[ p//2 ] >> ( (1-p%2)*4 ) ) & 15 ] )
                    bq_array.append( self.binaryqual[ p ] )
                    # now add any insertion later on
                    # get next cigar
                    if m + 1 == len( self.cigar ):
                        break
                    i = self.cigar[ m + 1 ]
                    op = i & 15
                    op_l = i >> 4
                    if op == 1:           #insertion
                        for n in range( op_l ):
                            p += 1
                            seq_array.append( __BAMDNACODE__[ (self.binaryseq[ p//2 ] >> ( (1-p%2)*4 ) ) & 15 ] )
                            bq_array.append( self.binaryqual[ p ] )
                        #print self.SEQ, seq_array
                    break
                else:
                    # go to the next cigar code
                    p += op_l
                    res -= op_l
            elif op in [ 2, 3 ]: # D N
                if res < op_l:
                    # find the position, however ...
                    # position located in a region in reference that not exists in query
                    pos = p
                    seq_array.append( b'*' )
                    bq_array.append( 93 )   #assign 93 for deletion
                    break
                else:
                    # go to the next cigar code
                    res -= op_l
            elif op == 1 :      # Insertion
                p += op_l
                # if res == 0:    # no residue left, so return a chunk of inserted sequence
                #     print "shouldn't run this code"
                #     # first, add the insertion point
                #     seq_array = bytearray( b'~' )
                #     bq_array.append( self.binaryqual[ p ] )
                #     # then add the inserted seq
                #     for i in range( op_l ):
                #         p += 1
                #         seq_array.append( __BAMDNACODE__[ (self.binaryseq[ p//2 ] >> ( (1-p%2)*4 ) ) & 15 ]  ) 
                #         bq_array.append( self.binaryqual[ p ] )
                #     break
                # else:
                #     p += op_l
            elif op == 4 :      # Softclip. If it's Softclip, we'd better not return the extra seq
                p += op_l

        if pos == 0 or pos == self.l - 1:
            tip = True
        else:
            tip = False
                
        return ( seq_array, bq_array, self.strand, tip, pos )
        # last position ?
        #raise Exception("Not expected to see this")

    cdef int relative_ref_pos_to_relative_query_pos ( self, long relative_ref_pos ):
        """Convert relative pos on ref to pos on query.
        """
        cdef:
            int p, res, op, op_l
        p = 0
        res = relative_ref_pos
        
        for i in self.cigar:
            op = i & 15
            op_l = i >> 4
            if op in [0, 7, 8]:         # M = X alignment match (match or mismatch)
                if res < op_l:
                    p += res
                    return p
                else:
                    p += op_l
                    res -= op_l
            elif op in [ 2, 3 ]: # D N
                if res < op_l:
                    # position located in a region in reference that not exists in query
                    return -1
                else:
                    res -= op_l
            elif op in [ 1, 4 ]:       # I
                p += op_l
        return p


### End ###

