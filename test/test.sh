python test_callsummits.py
macs2 callpeak -t random_test.bed -g 10000 -n random_test_CS --nomodel --extsize 150 --call-summits -B
macs2 callpeak -t random_test.bed -g 10000 -n random_test --nomodel --extsize 150
macs2 refinepeak -b random_test_CS_summits.bed -i random_test.bed -o random_test
bdg2bw random_test_CS_treat_pileup.bdg hg19.len
bedtools sort -i random_test.bed | bedtools bedtobam -i - -g hg19.len > random_test.bam
samtools index random_test.bam


# test all sub-commands

CHIP=ChIP_200K.bed.gz
CTRL=Control_200K.bed.gz

# callpeak

macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_narrow -B --outdir run_callpeak_narrow
macs2 callpeak -t $CHIP -c $CTRL -n run_callpeak_broad -B --outdir run_callpeak_broad --broad

# bdgcmp

macs2 bdgcmp -t run_callpeak_narrow/run_callpeak_narrow_treat_pileup.bdg -c run_callpeak_narrow/run_callpeak_narrow_control_lambda.bdg -m ppois FE --outdir run_bdgcmp_oprefix --o-prefix run_bdgcmp_oprefix

macs2 bdgcmp -t run_callpeak_narrow/run_callpeak_narrow_treat_pileup.bdg -c run_callpeak_narrow/run_callpeak_narrow_control_lambda.bdg -m ppois FE --outdir run_bdgcmp_ofile -o run_bdgcmp_ofile_ppois.bdg run_bdgcmp_ofile_FE.bdg

# bdgpeakcall

macs2 bdgpeakcall -i run_bdgcmp_single/run_bdgcmp_single_myFE.bdg -c 2 --outdir run_bdgpeakcall --o-prefix run_bdgpeakcall_w_prefix

macs2 bdgpeakcall -i run_bdgcmp_single/run_bdgcmp_single_myFE.bdg -c 2 --outdir run_bdgpeakcall -o run_bdgpeakcall_w_ofile.narrowPeak

# bdgbroadcall

macs2 bdgbroadcall -i run_bdgcmp_single/run_bdgcmp_single_myFE.bdg -c 4 -C 2 --outdir run_bdgbroadcall --o-prefix run_bdgbroadcall_w_prefix

macs2 bdgbroadcall -i run_bdgcmp_single/run_bdgcmp_single_myFE.bdg -c 4 -C 2 --outdir run_bdgbroadcall -o run_bdgbroadcall_w_ofile.gappedPeak

# diffpeak

#NA

# bdgdiff

macs2 callpeak -c $CHIP -t $CTRL -n run_callpeak_narrow_revert -B --outdir run_callpeak_narrow_revert

macs2 bdgdiff --t1 run_callpeak_narrow/run_callpeak_narrow_treat_pileup.bdg --c1 run_callpeak_narrow/run_callpeak_narrow_control_lambda.bdg --t2 run_callpeak_narrow_revert/run_callpeak_narrow_revert_treat_pileup.bdg --c2 run_callpeak_narrow_revert/run_callpeak_narrow_revert_control_lambda.bdg --o-prefix run_bdgdiff_prefix --outdir run_bdgdiff

macs2 bdgdiff --t1 run_callpeak_narrow/run_callpeak_narrow_treat_pileup.bdg --c1 run_callpeak_narrow/run_callpeak_narrow_control_lambda.bdg --t2 run_callpeak_narrow_revert/run_callpeak_narrow_revert_treat_pileup.bdg --c2 run_callpeak_narrow_revert/run_callpeak_narrow_revert_control_lambda.bdg -o cond1.bed cond2.bed common.bed --outdir run_bdgdiff

# filterdup

macs2 filterdup -i $CHIP --outdir run_filterdup -o run_filterdup_result.bed

# predictd

macs2 predictd -i $CHIP --outdir run_predictd --rfile run_predictd.R

# pileup

macs2 pileup -i $CHIP --outdir run_pileup -o run_pileup.bdg

# randsample

macs2 randsample -t $CHIP -n 100000 --seed 31415926 --outdir run_randsample -o run_randsample.bed

# refinepeak

macs2 refinepeak -b run_callpeak_narrow/run_callpeak_narrow_peaks.narrowPeak -i $CHIP --outdir run_refinepeak --o-prefix run_refinepeak_w_prefix

macs2 refinepeak -b run_callpeak_narrow/run_callpeak_narrow_peaks.narrowPeak -i $CHIP --outdir run_refinepeak -o run_refinepeak_w_ofile.bed

