# test all sub-commands

CHIP=CTCF_ChIP_200K.bed.gz
CTRL=CTCF_Control_200K.bed.gz

CHIPPE=CTCF_PE_ChIP_chr22.bam
CTRLPE=CTCF_PE_CTRL_chr22.bam

# callpeak
echo "callpeak"

macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow -B --outdir run_callpeak_narrow --call-summits &> run_callpeak_narrow/run_callpeak_narrow.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow2 -B --outdir run_callpeak_narrow --nomodel --extsize 100 &> run_callpeak_narrow/run_callpeak_narrow2.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow3 -B --outdir run_callpeak_narrow --nomodel --extsize 100 --shift -50 &> run_callpeak_narrow/run_callpeak_narrow3.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow4 -B --outdir run_callpeak_narrow --nomodel --nolambda --extsize 100 --shift -50 &> run_callpeak_narrow/run_callpeak_narrow4.log
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_broad -B --outdir run_callpeak_broad --broad &> run_callpeak_broad/run_callpeak_broad.log

echo "callpeak on PE"

macs2 callpeak -f BAMPE -t $CHIPPE -c $CTRLPE -n run_callpeak_pe_narrow -B --outdir run_callpeak_pe_narrow --call-summits &> run_callpeak_pe_narrow/run_callpeak_pe_narrow.log
macs2 callpeak -f BAMPE -t $CHIPPE -c $CTRLPE -n run_callpeak_pe_broad -B --outdir run_callpeak_pe_broad --broad &> run_callpeak_pe_broad/run_callpeak_pe_broad.log

# bdgcmp
echo "bdgcmp"

macs2 bdgcmp -t run_callpeak_narrow/run_callpeak_narrow_treat_pileup.bdg -c run_callpeak_narrow/run_callpeak_narrow_control_lambda.bdg -m ppois FE --outdir run_bdgcmp_oprefix --o-prefix run_bdgcmp_oprefix &> run_bdgcmp_oprefix/run_bdgcmp_oprefix.log

macs2 bdgcmp -t run_callpeak_narrow/run_callpeak_narrow_treat_pileup.bdg -c run_callpeak_narrow/run_callpeak_narrow_control_lambda.bdg -m ppois FE --outdir run_bdgcmp_ofile -o run_bdgcmp_ofile_ppois.bdg run_bdgcmp_ofile_FE.bdg &> run_bdgcmp_ofile/run_bdgcmp_ofile.log

# bdgpeakcall
echo "bdgpeakcall"

macs2 bdgpeakcall -i run_bdgcmp_single/run_bdgcmp_single_myFE.bdg -c 2 --outdir run_bdgpeakcall --o-prefix run_bdgpeakcall_w_prefix &> run_bdgpeakcall/run_bdgpeakcall_w_prefix.log

macs2 bdgpeakcall -i run_bdgcmp_single/run_bdgcmp_single_myFE.bdg -c 2 --outdir run_bdgpeakcall -o run_bdgpeakcall_w_ofile.narrowPeak &> run_bdgpeakcall/run_bdgpeakcall_w_ofile.log

# bdgbroadcall
echo "bdgbroadcall"

macs2 bdgbroadcall -i run_bdgcmp_single/run_bdgcmp_single_myFE.bdg -c 4 -C 2 --outdir run_bdgbroadcall --o-prefix run_bdgbroadcall_w_prefix &> run_bdgbroadcall/run_bdgbroadcall_w_prefix.log

macs2 bdgbroadcall -i run_bdgcmp_single/run_bdgcmp_single_myFE.bdg -c 4 -C 2 --outdir run_bdgbroadcall -o run_bdgbroadcall_w_ofile.gappedPeak &> run_bdgbroadcall/run_bdgbroadcall_w_ofile.log

# diffpeak

#NA

# bdgdiff
echo "bdgdiff"

macs2 callpeak -c $CHIP -t $CTRL -n run_callpeak_narrow_revert -B --outdir run_callpeak_narrow_revert &> run_callpeak_narrow_revert/run_callpeak_narrow_revert.log

macs2 bdgdiff --t1 run_callpeak_narrow/run_callpeak_narrow_treat_pileup.bdg --c1 run_callpeak_narrow/run_callpeak_narrow_control_lambda.bdg --t2 run_callpeak_narrow_revert/run_callpeak_narrow_revert_treat_pileup.bdg --c2 run_callpeak_narrow_revert/run_callpeak_narrow_revert_control_lambda.bdg --o-prefix run_bdgdiff_prefix --outdir run_bdgdiff &> run_bdgdiff/run_bdgdiff_w_prefix.log

macs2 bdgdiff --t1 run_callpeak_narrow/run_callpeak_narrow_treat_pileup.bdg --c1 run_callpeak_narrow/run_callpeak_narrow_control_lambda.bdg --t2 run_callpeak_narrow_revert/run_callpeak_narrow_revert_treat_pileup.bdg --c2 run_callpeak_narrow_revert/run_callpeak_narrow_revert_control_lambda.bdg -o cond1.bed cond2.bed common.bed --outdir run_bdgdiff &> run_bdgdiff/run_bdgdiff_w_o_file.log

# filterdup
echo "filterdup"

macs2 filterdup -i $CHIP --outdir run_filterdup -o run_filterdup_result.bed --dry-run &> run_filterdup/run_filterdup_d.log

macs2 filterdup -i $CHIP --outdir run_filterdup -o run_filterdup_result.bed &> run_filterdup/run_filterdup.log

# predictd
echo "predictd"

macs2 predictd -i $CHIP --outdir run_predictd --rfile run_predictd.R &> run_predictd/run_predictd.log

# pileup
echo "pileup"

macs2 pileup -i $CHIP --outdir run_pileup -o run_pileup.bdg &> run_pileup/run_pileup.log

# randsample
echo "randsample"

macs2 randsample -t $CHIP -n 100000 --seed 31415926 --outdir run_randsample -o run_randsample.bed &> run_randsample/run_randsample.log

# refinepeak
echo "refinepeak"

macs2 refinepeak -b run_callpeak_narrow/run_callpeak_narrow_peaks.narrowPeak -i $CHIP --outdir run_refinepeak --o-prefix run_refinepeak_w_prefix &> run_refinepeak/run_refinepeak_w_prefix.log

macs2 refinepeak -b run_callpeak_narrow/run_callpeak_narrow_peaks.narrowPeak -i $CHIP --outdir run_refinepeak -o run_refinepeak_w_ofile.bed &> run_refinepeak/run_refinepeak_w_ofile.log

