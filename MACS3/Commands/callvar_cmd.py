# Time-stamp: <2022-09-15 17:25:49 Tao Liu>

"""Description: Call variants directly

Copyright (c) 2017-2023 Tao Liu <vladimir.liu@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status: release candidate
@version: $Id$
@author:  Tao Liu
@contact: vladimir.liu@gmail.com
"""

# ------------------------------------
# python modules
# ------------------------------------

import datetime
import sys
from functools import partial
import multiprocessing as mp

from time import time
from math import ceil

# ------------------------------------
# own python modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Utilities.OptValidator import opt_validate_callvar
from MACS3.IO.PeakIO import PeakIO
from MACS3.IO.BAM import BAMaccessor
from MACS3.Signal.RACollection import RACollection
from MACS3.Signal.PeakVariants import PeakVariants


VCFHEADER_0="""##fileformat=VCFv4.1
##fileDate=%s
##source=MACS_V%s
##Program_Args=%s
##INFO=<ID=M,Number=.,Type=String,Description="MACS Model with minimum BIC value">
##INFO=<ID=MT,Number=.,Type=String,Description="Mutation type: SNV/Insertion/Deletion">
##INFO=<ID=DPT,Number=1,Type=Integer,Description="Depth Treatment: Read depth in ChIP-seq data">
##INFO=<ID=DPC,Number=1,Type=Integer,Description="Depth Control: Read depth in control data">
##INFO=<ID=DP1T,Number=.,Type=String,Description="Read depth of top1 allele in ChIP-seq data">
##INFO=<ID=DP2T,Number=.,Type=String,Description="Read depth of top2 allele in ChIP-seq data">
##INFO=<ID=DP1C,Number=.,Type=String,Description="Read depth of top1 allele in control data">
##INFO=<ID=DP2C,Number=.,Type=String,Description="Read depth of top2 allele in control data">
##INFO=<ID=lnLHOMOMAJOR,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of homozygous with major allele model">
##INFO=<ID=lnLHOMOMINOR,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of homozygous with minor allele model">
##INFO=<ID=lnLHETERNOAS,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of heterozygous with no allele-specific model">
##INFO=<ID=lnLHETERAS,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of heterozygous with allele-specific model">
##INFO=<ID=BICHOMOMAJOR,Number=1,Type=Float,Description="BIC value of homozygous with major allele model">
##INFO=<ID=BICHOMOMINOR,Number=1,Type=Float,Description="BIC value of homozygous with minor allele model">
##INFO=<ID=BICHETERNOAS,Number=1,Type=Float,Description="BIC value of heterozygous with no allele-specific model">
##INFO=<ID=BICHETERAS,Number=1,Type=Float,Description="BIC value of heterozygous with allele-specific model">
##INFO=<ID=GQHOMO,Number=1,Type=Integer,Description="Genotype quality of homozygous with major allele model">
##INFO=<ID=GQHETERNOAS,Number=1,Type=Integer,Description="Genotype quality of heterozygous with no allele-specific model">
##INFO=<ID=GQHETERAS,Number=1,Type=Integer,Description="Genotype quality of heterozygous with allele-specific model">
##INFO=<ID=GQHETERASsig,Number=1,Type=Integer,Description="Genotype quality of allele-specific significance compared with no allele-specific model">
##INFO=<ID=AR,Number=1,Type=Float,Description="Estimated allele ratio of heterozygous with allele-specific model">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="read depth after filtering">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality score">"""

VCFHEADER="""##fileformat=VCFv4.1
##fileDate=%s
##source=MACS_V%s
##Program_Args=%s
##INFO=<ID=M,Number=.,Type=String,Description="MACS Model with minimum BIC value">
##INFO=<ID=MT,Number=.,Type=String,Description="Mutation type: SNV/Insertion/Deletion">
##INFO=<ID=DPT,Number=1,Type=Integer,Description="Depth Treatment: Read depth in ChIP-seq data">
##INFO=<ID=DPC,Number=1,Type=Integer,Description="Depth Control: Read depth in control data">
##INFO=<ID=DP1T,Number=.,Type=String,Description="Read depth of top1 allele in ChIP-seq data">
##INFO=<ID=DP2T,Number=.,Type=String,Description="Read depth of top2 allele in ChIP-seq data">
##INFO=<ID=DP1C,Number=.,Type=String,Description="Read depth of top1 allele in control data">
##INFO=<ID=DP2C,Number=.,Type=String,Description="Read depth of top2 allele in control data">
##INFO=<ID=DBIC,Number=.,Type=Float,Description="Difference of BIC of selected model vs second best alternative model">
##INFO=<ID=BICHOMOMAJOR,Number=1,Type=Integer,Description="BIC of homozygous with major allele model">
##INFO=<ID=BICHOMOMINOR,Number=1,Type=Integer,Description="BIC of homozygous with minor allele model">
##INFO=<ID=BICHETERNOAS,Number=1,Type=Integer,Description="BIC of heterozygous with no allele-specific model">
##INFO=<ID=BICHETERAS,Number=1,Type=Integer,Description="BIC of heterozygous with allele-specific model">
##INFO=<ID=AR,Number=1,Type=Float,Description="Estimated allele ratio of heterozygous with allele-specific model">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth after filtering bad reads">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality score">
##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled genotype likelihoods for 00, 01, 11 genotype">"""

# ------------------------------------
# Main function
# ------------------------------------

def check_names(treat, control, error_stream):
    """check common chromosome names"""
    tchrnames = set(treat.get_chr_names())
    cchrnames = set(control.get_chr_names())
    commonnames = tchrnames.intersection(cchrnames)
    if len(commonnames)==0:
        error_stream("No common chromosome names can be found from treatment and control! Check your input files! MACS will quit...")
        error_stream("Chromosome names in treatment: %s" % ",".join(sorted(tchrnames)))
        error_stream("Chromosome names in control: %s" % ",".join(sorted(cchrnames)))
        sys.exit()


def run( args ):
    """The Main function/pipeline for MACS

    """
    options = opt_validate_callvar( args )

    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    peakbedfile = options.peakbed
    tfile = options.tfile
    cfile = options.cfile
    top2allelesminr = options.top2allelesMinRatio
    min_altallele_count = options.altalleleMinCount
    max_allowed_ar = options.maxAR
    NP = options.np
    if NP<=0:
        NP = 1
    min_homo_GQ = options.GQCutoffHomo
    min_heter_GQ = options.GQCutoffHetero
    minQ = options.Q
    maxDuplicate = options.maxDuplicate

    # parameter for assembly
    fermiMinOverlap = options.fermiMinOverlap
    fermi = options.fermi
    
    peakio = open( peakbedfile )
    peaks = PeakIO()
    i = 0
    for l in peakio:
        fs = l.rstrip().split()
        i += 1
        peaks.add( fs[0].encode(), int(fs[1]), int(fs[2]), name=b"%d" % i )
    peaks.sort()

    chrs = peaks.get_chr_names()

    tbam = BAMaccessor( tfile )
    if cfile:
        cbam = BAMaccessor( cfile )
        assert tbam.get_chromosomes()[0] in cbam.get_chromosomes() or cbam.get_chromosomes()[0] in tbam.get_chromosomes(), Exception("It seems Treatment and Control BAM use different naming for chromosomes! Check headers of both files.")
        #assert tbam.get_chromosomes() == cbam.get_chromosomes(), Exception("Treatment and Control BAM files have different orders of sorted chromosomes! Please check BAM Headers and re-sort BAM files.")
    else:
        cbam = None


    ra_collections = []

    # prepare and write header of output file (.vcf)
    ovcf = open(options.ofile, "w")
    tmpcmdstr = " --fermi "+ fermi+ " --fermi-overlap "+str(fermiMinOverlap)
    ovcf.write ( VCFHEADER % (datetime.date.today().strftime("%Y%m%d"), MACS_VERSION, " ".join(sys.argv[1:] + ["-Q", str(minQ), "-D", str(maxDuplicate), "--max-ar", str(max_allowed_ar), "--top2alleles-mratio", str(top2allelesminr), "--top2allele-count", str(min_altallele_count), "-g", str(min_heter_GQ), "-G", str(min_homo_GQ), tmpcmdstr]) ) + "\n" )
    for (chrom, chrlength) in tbam.get_rlengths().items():
        ovcf.write( "##contig=<ID=%s,length=%d,assembly=NA>\n" % ( chrom.decode(), chrlength ) )
    ovcf.write ( "\t".join( ("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE") ) + "\n" )

    # to get time
    t_total = 0
    t_prepare_ra = 0
    t_assemble = 0
    #t_call_top2alleles = 0
    #t_call_lnL = 0
    t_call_variants = 0
    t_call_GT = 0
    #t_call_to_vcf = 0
    t_total_0 = time()

    for chrom in sorted(tbam.get_chromosomes()):
        peaks_chr = peaks.get_data_from_chrom( chrom )
        for peak in peaks_chr:
            # note, when we extract reads from BAM within a peak
            # region, we assume BAM should be sorted and the BAM
            # should be generated from "samtools view -L" process.

            # print ( "---begin of peak---")
            info ( f"Peak: {chrom.decode()} {peak['start']} {peak['end']}" )

            flag_todo_lassembly = False

            t0 = time()
            try:
                if cbam:
                    ra_collection = RACollection( chrom, peak, tbam.get_reads_in_region( chrom, peak["start"], peak["end"], maxDuplicate=maxDuplicate ), cbam.get_reads_in_region( chrom, peak["start"], peak["end"], maxDuplicate=maxDuplicate) )
                else:
                    ra_collection = RACollection( chrom, peak, tbam.get_reads_in_region( chrom, peak["start"], peak["end"], maxDuplicate=maxDuplicate ) )
            except:
                info ("No reads found in peak: ", chrom.decode(), peak["start"], peak["end"], ". Skipped!")
                # while there is no reads in peak region, simply skip it.
                continue
            
            ra_collection.remove_outliers( percent = 5 )
            t_prepare_ra += time() - t0

            # print ( "Reads in Peak:")
            # print ( ra_collection.get_FASTQ().decode() )

            s = ra_collection["peak_refseq"]

            peak_variants = PeakVariants( chrom.decode(), peak["start"], peak["end"], s )
            

            if fermi == "auto" or fermi == "off":
                # first pass to call variant w/o assembly
                # multiprocessing the following part
                t_call_variants_0 = time()
                info ( " Call variants w/o assembly")

                # -- now make multi processes
                # divide right-left into NP parts
                window_size = ceil( ( ra_collection["right"] - ra_collection["left"] ) / NP )

                P = mp.Pool( NP )
                # this partial function will only be used in multiprocessing
                p_call_variants_at_range =  partial(call_variants_at_range, s=s, collection=ra_collection, top2allelesminr=top2allelesminr, max_allowed_ar = max_allowed_ar, min_altallele_count = min_altallele_count, min_homo_GQ = min_homo_GQ, min_heter_GQ = min_heter_GQ, minQ=minQ)

                ranges = []
                for i in range( NP ):
                    l = i * window_size + ra_collection["left"]
                    r = min( (i + 1) * window_size + ra_collection["left"], ra_collection["right"] )
                    ranges.append( (l, r) )

                mapresults = P.map_async( p_call_variants_at_range, ranges )
                P.close()
                P.join()
                results = mapresults.get(timeout=window_size*300)
                for i in range( NP ):
                    for result in results[ i ]:
                        peak_variants.add_variant( result[0], result[1] )

                t_call_variants += time() - t_call_variants_0

            # Next, check if we should do local assembly
            if ( fermi == "auto" and ( peak_variants.has_indel() or peak_variants.has_refer_biased_01() ) ) or fermi == "on":
                #print( peak_variants.has_indel() )
                #print( peak_variants.has_refer_biased_01() )
                    
                # invoke fermi to assemble local sequence and filter out those can not be mapped to unitigs.
                info ( " Try to call variants w/ fermi-lite assembly")
                unitig_collection = ra_collection.build_unitig_collection( fermiMinOverlap )
                if unitig_collection == -1:
                    info(" Too many mismatches found while assembling the sequence, we will skip this region entirely!")
                    continue
                elif unitig_collection == 0:
                    info ( "  Failed to assemble unitigs, fall back to previous results" )
                    if peak_variants.n_variants() > 0:
                        peak_variants.fix_indels()
                        ovcf.write( peak_variants.toVCF() )
                        continue
                # uncomment the following to print those assembled unitigs and their alignments to reference genome
                #for u in unitig_collection["URAs_list"]:
                #    print( u["seq"].decode(), u["lpos"], u["rpos"], u["count"] )
                #    print( "a",u["unitig_aln"].decode() )
                #    print( "r",u["reference_aln"].decode() )
            else:
                # if we do not assemble, write results now
                if peak_variants.n_variants() > 0:
                    peak_variants.fix_indels()                    
                    ovcf.write( peak_variants.toVCF() )
                    continue

            # reach here only if we need l assembly and the assembly returns result

            # If a peak has no indel but with refer_biased_01, we
            # revisit all refer_biased_01 now. We do not use
            # multiprocessing here for simplicity since there won't be
            # too many in a peak region.
            if ( fermi == "auto" and ( not peak_variants.has_indel() ) and peak_variants.has_refer_biased_01()  ):
                pos_tobe_revisit = peak_variants.get_refer_biased_01s()
                for i in pos_tobe_revisit:
                    ref_nt = chr(s[ i-ra_collection["left"] ] ).encode()
                    if ref_nt == b'N':
                        peak_variants.remove_variant( i )
                        continue
                    PRI = unitig_collection.get_PosReadsInfo_ref_pos ( i, ref_nt, Q=minQ )
                    if PRI.raw_read_depth( opt="T" ) == 0: # skip if coverage is 0
                        peak_variants.remove_variant( i )                        
                        continue
                    PRI.update_top_alleles( top2allelesminr, min_altallele_count, max_allowed_ar )
                    PRI.call_GT( max_allowed_ar )
                    PRI.apply_GQ_cutoff(min_homo_GQ, min_heter_GQ)
                    if not PRI.filterflag():
                        peak_variants.replace_variant( i, PRI.toVariant() )
                    else:
                        peak_variants.remove_variant( i )
                if peak_variants.n_variants() > 0:
                    peak_variants.fix_indels()
                    ovcf.write( peak_variants.toVCF() )
                    continue

            # in this case, we call variants at every locations in the peak based on local assembly.
            if ( fermi == "on" or ( fermi == "auto" and peak_variants.has_indel() ) ):
                peak_variants = PeakVariants( chrom.decode(), peak["start"], peak["end"], s ) #reset

                # --- make multi processes
                # divide right-left into NP parts
                window_size = ceil( ( ra_collection["right"] - ra_collection["left"] ) / NP )
                P = mp.Pool( NP )
                # this partial function will only be used in multiprocessing
                p_call_variants_at_range =  partial(call_variants_at_range, s=s, collection=unitig_collection, top2allelesminr=top2allelesminr, max_allowed_ar = max_allowed_ar, min_altallele_count = min_altallele_count, min_homo_GQ = min_homo_GQ, min_heter_GQ = min_heter_GQ, minQ=minQ)
                ranges = []
                for i in range( NP ):
                    l = i * window_size + ra_collection["left"]
                    r = min( (i + 1) * window_size + ra_collection["left"], ra_collection["right"] )
                    ranges.append( (l, r) )

                mapresults = P.map_async( p_call_variants_at_range, ranges )
                P.close()
                P.join()
                results = mapresults.get(timeout=window_size*300)                
                for i in range( NP ):
                    for result in results[ i ]:
                        peak_variants.add_variant( result[0], result[1] )
                                
            if peak_variants.n_variants() > 0:
                peak_variants.fix_indels()
                ovcf.write( peak_variants.toVCF() )

    #print ("time to retrieve read alignment information from BAM:",t_prepare_ra,"(",round( 100 * t_prepare_ra/t_total, 2),"% )")
    return

def call_variants_at_range ( lr, s, collection, top2allelesminr, max_allowed_ar, min_altallele_count, min_homo_GQ, min_heter_GQ, minQ ):
#def call_variants_at_range ( lr, chrom, s, collection, top2allelesminr, max_allowed_ar, min_homo_GQ, min_heter_GQ ):
    result = []
    for i in range( lr[ 0 ], lr[ 1 ] ):
        ref_nt = chr(s[ i-collection["left"] ] ).encode()
        if ref_nt == b'N':
            continue

        PRI = collection.get_PosReadsInfo_ref_pos ( i, ref_nt, Q=minQ )
        if PRI.raw_read_depth( opt="T" ) == 0: # skip if coverage is 0
            continue

        PRI.update_top_alleles( top2allelesminr, min_altallele_count, max_allowed_ar )
        if not PRI.filterflag():
            #PRI.update_top_alleles( top2allelesminr )
            PRI.call_GT( max_allowed_ar )
            PRI.apply_GQ_cutoff(min_homo_GQ, min_heter_GQ)
        if not PRI.filterflag():
            result.append( ( i, PRI.toVariant() ) )
    return result


