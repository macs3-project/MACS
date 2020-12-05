#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "fml.h"

int main(int argc, char *argv[])
{
	fml_opt_t opt;
	int c, n_seqs, n_utg, gfa_out = 0;
	bseq1_t *seqs;
	fml_utg_t *utg;

	fml_opt_init(&opt);
	while ((c = getopt(argc, argv, "gOAe:l:r:t:c:d:v:")) >= 0) {
		if (c == 'e') opt.ec_k = atoi(optarg);
		else if (c == 'l') opt.min_asm_ovlp = atoi(optarg);
		else if (c == 'r') opt.mag_opt.min_dratio1 = atof(optarg);
		else if (c == 'A') opt.mag_opt.flag |= MAG_F_AGGRESSIVE;
		else if (c == 'O') opt.mag_opt.flag &= ~MAG_F_POPOPEN;
		else if (c == 'd') opt.mag_opt.max_bdiff = atoi(optarg);
		else if (c == 't') opt.n_threads = atoi(optarg);
		else if (c == 'g') gfa_out = 1;
		else if (c == 'v') fm_verbose = atoi(optarg);
		else if (c == 'c') {
			char *p;
			opt.min_cnt = strtol(optarg, &p, 10);
			if (*p == ',') opt.max_cnt = strtol(p + 1, &p, 10);
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fml-asm [options] <in.fq>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -e INT          k-mer length for error correction (0 for auto; -1 to disable) [%d]\n", opt.ec_k);
		fprintf(stderr, "  -c INT1[,INT2]  range of k-mer & read count thresholds for ec and graph cleaning [%d,%d]\n", opt.min_cnt, opt.max_cnt);
		fprintf(stderr, "  -l INT          min overlap length during initial assembly [%d]\n", opt.min_asm_ovlp);
		fprintf(stderr, "  -r FLOAT        drop an overlap if its length is below maxOvlpLen*FLOAT [%g]\n", opt.mag_opt.min_dratio1);
		fprintf(stderr, "  -t INT          number of threads (don't use multi-threading for small data sets) [%d]\n", opt.n_threads);
		fprintf(stderr, "  -d INT          retain a bubble if one side is longer than the other side by >INT-bp [%d]\n", opt.mag_opt.max_bdiff);
		fprintf(stderr, "  -A              discard heterozygotes (apply this to assemble bacterial genomes; override -O)\n");
		fprintf(stderr, "  -O              don't apply aggressive tip trimming\n");
		fprintf(stderr, "  -g              output the assembly graph in the GFA format\n");
		return 1;
	}
	seqs = bseq_read(argv[optind], &n_seqs);
	utg = fml_assemble(&opt, n_seqs, seqs, &n_utg);
	if (!gfa_out) fml_utg_print(n_utg, utg);
	else fml_utg_print_gfa(n_utg, utg);
	fml_utg_destroy(n_utg, utg);
	return 0;
}
