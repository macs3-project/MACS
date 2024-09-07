# VCF format

The VCF format is used as the output from `callvar` subcommand. The
version of format we use in MACS3 is
[v4.1](https://samtools.github.io/hts-specs/VCFv4.1.pdf). Please check
the link for detail. 

The output from `callvar` has the following customized fields as
defined in the VCF header lines such as:

```
##fileformat=VCFv4.1
##fileDate=20240514
##source=MACS_V3.0.1
##Program_Args=callvar -b callvar_testing.narrowPeak -t CTCF_PE_ChIP_chr22_50k.bam -c CTCF_PE_CTRL_chr22_50k.bam -o ../temp/test513_run_callvar/PEsample.vcf -Q 20 -D 1 --max-ar 0.95 --top2alleles-mratio 0.8 --top2allele-count 2 -g 0 -G 0  --fermi auto --fermi-overlap 30
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
##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled genotype likelihoods for 00, 01, 11 genotype"
```

The header lines contain the command line options used to generate
this output, date of the file, and definitions of the customized
fields in 'INFO' and 'FORMAT'/'SAMPLE' columns of the VCF. Here is an
example of the actual data:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr22	17255752	.	A	G	58	.	M=heter_unsure;MT=SNV;DPT=7;DPC=0;DP1T=5A;DP2T=2G;DP1C=0A;DP2C=0G;SB=0,0,5,2;DBIC=23.21;BICHOMOMAJOR=37.77;BICHOMOMINOR=84.27;BICHETERNOAS=13.53;BICHETERAS=14.56;AR=0.71	GT:DP:GQ:PL	0/1:7:58:159,0,58
chr22	17392539	.	G	C	138	.	M=heter_noAS;MT=SNV;DPT=13;DPC=0;DP1T=7C;DP2T=6G;DP1C=0C;DP2C=0G;SB=0,1,7,5;DBIC=61.11;BICHOMOMAJOR=84.75;BICHOMOMINOR=101.33;BICHETERNOAS=23.63;BICHETERAS=26.12;AR=0.54	GT:DP:GQ:PL	0/1:13:138:174,0,13
```

