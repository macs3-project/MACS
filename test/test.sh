#!/bin/bash

# integrative subcmds testing

if [ $# -lt 1 ];then
    echo "Run all tests for subcommands of MACS2. Need 1 parameter for a tag name! A unique string combining date, time and MACS2 version is recommended. ./test.sh <TAG>"
    exit
fi

# test all sub-commands
TAG=$1

CHIP=CTCF_ChIP_200K.bed.gz
CTRL=CTCF_Control_200K.bed.gz

CHIPPE=CTCF_PE_ChIP_chr22.bam
CTRLPE=CTCF_PE_CTRL_chr22.bam

CHIPBEDPE=CTCF_PE_ChIP_chr22.bedpe
CTRLBEDPE=CTCF_PE_CTRL_chr22.bedpe

# callpeak
echo "callpeak"

mkdir ${TAG}_run_callpeak_narrow

macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow0 -B --outdir ${TAG}_run_callpeak_narrow  &> ${TAG}_run_callpeak_narrow/run_callpeak_narrow0.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow1 -B --outdir ${TAG}_run_callpeak_narrow --d-min 15 --call-summits &> ${TAG}_run_callpeak_narrow/run_callpeak_narrow1.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow2 -B --outdir ${TAG}_run_callpeak_narrow --nomodel --extsize 100 &> ${TAG}_run_callpeak_narrow/run_callpeak_narrow2.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow3 -B --outdir ${TAG}_run_callpeak_narrow --nomodel --extsize 100 --shift -50 &> ${TAG}_run_callpeak_narrow/run_callpeak_narrow3.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow4 -B --outdir ${TAG}_run_callpeak_narrow --nomodel --nolambda --extsize 100 --shift -50 &> ${TAG}_run_callpeak_narrow/run_callpeak_narrow4.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow5 -B --outdir ${TAG}_run_callpeak_narrow --scale-to large &> ${TAG}_run_callpeak_narrow/run_callpeak_narrow5.log

mkdir ${TAG}_run_callpeak_broad 

macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_broad -B --outdir ${TAG}_run_callpeak_broad --broad &> ${TAG}_run_callpeak_broad/run_callpeak_broad.log

echo "callpeak on PE"

mkdir  ${TAG}_run_callpeak_pe_narrow
mkdir  ${TAG}_run_callpeak_pe_broad

macs2 callpeak -f BAMPE -t $CHIPPE -c $CTRLPE -n run_callpeak_bampe_narrow -B --outdir ${TAG}_run_callpeak_pe_narrow --call-summits &> ${TAG}_run_callpeak_pe_narrow/run_callpeak_bampe_narrow.log
macs2 callpeak -f BAMPE -t $CHIPPE -c $CTRLPE -n run_callpeak_bampe_broad -B --outdir ${TAG}_run_callpeak_pe_broad --broad &> ${TAG}_run_callpeak_pe_broad/run_callpeak_bampe_broad.log
macs2 callpeak -f BEDPE -t $CHIPBEDPE -c $CTRLBEDPE -n run_callpeak_bedpe_narrow -B --outdir ${TAG}_run_callpeak_pe_narrow --call-summits &> ${TAG}_run_callpeak_pe_narrow/run_callpeak_bedpe_narrow.log
macs2 callpeak -f BEDPE -t $CHIPBEDPE -c $CTRLBEDPE -n run_callpeak_bedpe_broad -B --outdir ${TAG}_run_callpeak_pe_broad --broad &> ${TAG}_run_callpeak_pe_broad/run_callpeak_bedpe_broad.log
macs2 callpeak -f BEDPE -t $CHIPBEDPE -n run_callpeak_pe_narrow_onlychip -B --outdir ${TAG}_run_callpeak_pe_narrow &> ${TAG}_run_callpeak_pe_narrow/run_callpeak_pe_narrow_onlychip.log

# pileup
echo "pileup"

mkdir ${TAG}_run_pileup

macs2 pileup -f BED -i $CHIP --extsize 200 --outdir ${TAG}_run_pileup -o run_pileup_ChIP.bed.bdg &> ${TAG}_run_pileup/run_pileup_ChIP.bed.log
macs2 pileup -f BED -i $CTRL --extsize 200 --outdir ${TAG}_run_pileup -o run_pileup_CTRL.bed.bdg &> ${TAG}_run_pileup/run_pileup_CTRL.bed.log

macs2 pileup -f BAMPE -i $CHIPPE --outdir ${TAG}_run_pileup -o run_pileup_ChIPPE.bampe.bdg &> ${TAG}_run_pileup/run_pileup_ChIPPE.bampe.log
macs2 pileup -f BAMPE -i $CTRLPE --outdir ${TAG}_run_pileup -o run_pileup_CTRLPE.bampe.bdg &> ${TAG}_run_pileup/run_pileup_CTRLPE.bampe.log

macs2 pileup -f BEDPE -i $CHIPBEDPE --outdir ${TAG}_run_pileup -o run_pileup_ChIPPE.bedpe.bdg &> ${TAG}_run_pileup/run_pileup_ChIPPE.bedpe.log
macs2 pileup -f BEDPE -i $CTRLBEDPE --outdir ${TAG}_run_pileup -o run_pileup_CTRLPE.bedpe.bdg &> ${TAG}_run_pileup/run_pileup_CTRLPE.bedpe.log

# filterdup
echo "filterdup"

mkdir ${TAG}_run_filterdup 

macs2 filterdup -i $CHIP --outdir ${TAG}_run_filterdup -o run_filterdup_result.bed --dry-run &> ${TAG}_run_filterdup/run_filterdup_d.log
macs2 filterdup -i $CHIP --outdir ${TAG}_run_filterdup -o run_filterdup_result.bed &> ${TAG}_run_filterdup/run_filterdup.log

macs2 filterdup -i $CHIPPE -f BAMPE --outdir ${TAG}_run_filterdup -o run_filterdup_result_pe.bedpe --dry-run &> ${TAG}_run_filterdup/run_filterdup_pe_d.log
macs2 filterdup -i $CHIPPE -f BAMPE --outdir ${TAG}_run_filterdup -o run_filterdup_result_pe.bedpe &> ${TAG}_run_filterdup/run_filterdup_pe.log

# predictd
echo "predictd"

mkdir ${TAG}_run_predictd

macs2 predictd -i $CHIP --d-min 10 --outdir ${TAG}_run_predictd --rfile run_predictd.R &> ${TAG}_run_predictd/run_predictd.log

# randsample
echo "randsample"

mkdir ${TAG}_run_randsample

macs2 randsample -i $CHIP -n 100000 --seed 31415926 --outdir ${TAG}_run_randsample -o run_randsample.bed &> ${TAG}_run_randsample/run_randsample.log

macs2 randsample -f BAMPE -i $CHIPPE -n 100000 --seed 31415926 --outdir ${TAG}_run_randsample -o run_randsample_bampe.bedpe &> ${TAG}_run_randsample/run_randsample_bampe.log

macs2 randsample -f BEDPE -i $CHIPBEDPE -n 100000 --seed 31415926 --outdir ${TAG}_run_randsample -o run_randsample_bedpe.bedpe &> ${TAG}_run_randsample/run_randsample_bedpe.log

# refinepeak
echo "refinepeak"

mkdir ${TAG}_run_refinepeak

macs2 refinepeak -b ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_peaks.narrowPeak -i $CHIP --outdir ${TAG}_run_refinepeak --o-prefix run_refinepeak_w_prefix &> ${TAG}_run_refinepeak/run_refinepeak_w_prefix.log

macs2 refinepeak -b ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_peaks.narrowPeak -i $CHIP --outdir ${TAG}_run_refinepeak -o run_refinepeak_w_ofile.bed &> ${TAG}_run_refinepeak/run_refinepeak_w_ofile.log

# bdgcmp
echo "bdgcmp"

mkdir ${TAG}_run_bdgcmp

macs2 bdgcmp -t ${TAG}_run_pileup/run_pileup_ChIP.bed.bdg -c ${TAG}_run_pileup/run_pileup_CTRL.bed.bdg  -m ppois FE -p 1 --outdir ${TAG}_run_bdgcmp --o-prefix run_bdgcmp &> ${TAG}_run_bdgcmp/run_bdgcmp.log

# bdgpeakcall
echo "bdgpeakcall"

mkdir ${TAG}_run_bdgpeakcall

macs2 bdgpeakcall -i ${TAG}_run_bdgcmp/run_bdgcmp_FE.bdg -c 2 --outdir ${TAG}_run_bdgpeakcall --o-prefix run_bdgpeakcall_w_prefix &> ${TAG}_run_bdgpeakcall/run_bdgpeakcall_w_prefix.log

# bdgbroadcall
echo "bdgbroadcall"

mkdir ${TAG}_run_bdgbroadcall

macs2 bdgbroadcall -i ${TAG}_run_bdgcmp/run_bdgcmp_FE.bdg -c 2 -C 1.5 --outdir ${TAG}_run_bdgbroadcall --o-prefix run_bdgbroadcall_w_prefix &> ${TAG}_run_bdgbroadcall/run_bdgbroadcall_w_prefix.log

# bdgdiff
echo "bdgdiff"

mkdir ${TAG}_run_callpeak_narrow_revert
mkdir ${TAG}_run_bdgdiff

macs2 callpeak -c $CHIP -t $CTRL -n run_callpeak_narrow_revert -B --outdir ${TAG}_run_callpeak_narrow_revert &> ${TAG}_run_callpeak_narrow_revert/run_callpeak_narrow_revert.log

macs2 bdgdiff --t1 ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_treat_pileup.bdg --c1 ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_control_lambda.bdg --t2 ${TAG}_run_callpeak_narrow_revert/run_callpeak_narrow_revert_treat_pileup.bdg --c2 ${TAG}_run_callpeak_narrow_revert/run_callpeak_narrow_revert_control_lambda.bdg --o-prefix run_bdgdiff_prefix --outdir ${TAG}_run_bdgdiff &> ${TAG}_run_bdgdiff/run_bdgdiff_w_prefix.log

macs2 bdgdiff --t1 ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_treat_pileup.bdg --c1 ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_control_lambda.bdg --t2 ${TAG}_run_callpeak_narrow_revert/run_callpeak_narrow_revert_treat_pileup.bdg --c2 ${TAG}_run_callpeak_narrow_revert/run_callpeak_narrow_revert_control_lambda.bdg -o cond1.bed cond2.bed common.bed --outdir ${TAG}_run_bdgdiff &> ${TAG}_run_bdgdiff/run_bdgdiff_w_o_file.log

# cmbreps
echo "cmbreps"

mkdir ${TAG}_run_cmbreps

macs2 cmbreps -i ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_treat_pileup.bdg ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_control_lambda.bdg ${TAG}_run_bdgcmp/run_bdgcmp_ppois.bdg -m max -o ${TAG}_cmbreps_max.bdg --outdir ${TAG}_run_cmbreps &> ${TAG}_run_cmbreps/run_cmbreps_max.log
macs2 cmbreps -i ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_treat_pileup.bdg ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_control_lambda.bdg ${TAG}_run_bdgcmp/run_bdgcmp_ppois.bdg -m mean -o ${TAG}_cmbreps_mean.bdg --outdir ${TAG}_run_cmbreps &> ${TAG}_run_cmbreps/run_cmbreps_mean.log
macs2 cmbreps -i ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_treat_pileup.bdg ${TAG}_run_callpeak_narrow/run_callpeak_narrow0_control_lambda.bdg ${TAG}_run_bdgcmp/run_bdgcmp_ppois.bdg -m fisher -o ${TAG}_cmbreps_fisher.bdg --outdir ${TAG}_run_cmbreps &> ${TAG}_run_cmbreps/run_cmbreps_fisher.log
