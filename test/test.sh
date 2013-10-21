python test_callsummits.py
macs2 callpeak -t random_test.bed -g 10000 -n random_test_CS --nomodel --extsize 150 --call-summits -B
macs2 callpeak -t random_test.bed -g 10000 -n random_test --nomodel --extsize 150
macs2 refinepeak -b random_test_CS_summits.bed -i random_test.bed -o random_test
bdg2bw random_test_CS_treat_pileup.bdg hg19.len
bedtools sort -i random_test.bed | bedtools bedtobam -i - -g hg19.len > random_test.bam
samtools index random_test.bam
