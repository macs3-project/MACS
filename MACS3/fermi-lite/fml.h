#ifndef FML_H
#define FML_H

#define FML_VERSION "r53"

#include <stdint.h>

typedef struct {
	int32_t l_seq;
	char *seq, *qual; // NULL-terminated strings; length expected to match $l_seq
} bseq1_t;

#define MAG_F_AGGRESSIVE 0x20 // pop variant bubbles (not default)
#define MAG_F_POPOPEN    0x40 // aggressive tip trimming (default)
#define MAG_F_NO_SIMPL   0x80 // skip bubble simplification (default)

typedef struct {
	int flag, min_ovlp, min_elen, min_ensr, min_insr, max_bdist, max_bdiff, max_bvtx, min_merge_len, trim_len, trim_depth;
	float min_dratio1, max_bcov, max_bfrac;
} magopt_t;

typedef struct {
	int n_threads;        // number of threads; don't use multi-threading for small data sets
	int ec_k;             // k-mer length for error correction; 0 for auto estimate
	int min_cnt, max_cnt; // both occ threshold in ec and tip threshold in cleaning lie in [min_cnt,max_cnt]
	int min_asm_ovlp;     // min overlap length during assembly
	int min_merge_len;    // during assembly, don't explicitly merge an overlap if shorter than this value
	magopt_t mag_opt;     // graph cleaning options
} fml_opt_t;

struct rld_t;
struct mag_t;

typedef struct {
        // uint32_t len_:31, from_:1; // $from and $to: 0 meaning overlapping 5'-end; 1 overlapping 3'-end
        // uint32_t id_:31, to_:1;    // $id: unitig number
        uint32_t len_;
        uint32_t from_;
        uint32_t id_;
        uint32_t to_;
} fml_ovlp_t;

typedef struct {
	int32_t len;      // length of sequence
	int32_t nsr;      // number of supporting reads
	char *seq;        // unitig sequence
	char *cov;        // cov[i]-33 gives per-base coverage at i
	int n_ovlp[2];    // number of 5'-end [0] and 3'-end [1] overlaps
	fml_ovlp_t *ovlp; // overlaps, of size n_ovlp[0]+n_ovlp[1]
} fml_utg_t;

extern int fm_verbose;

#ifdef __cplusplus
extern "C" {
#endif

/************************
 * High-level functions *
 ************************/

/**
 * Read all sequences from a FASTA/FASTQ file
 *
 * @param fn       filename; NULL or "-" for stdin
 * @param n        (out) number of sequences read into RAM
 *
 * @return array of sequences
 */
bseq1_t *bseq_read(const char *fn, int *n);

/**
 * Initialize default parameters
 *
 * @param opt      (out) pointer to parameters
 */
void fml_opt_init(fml_opt_t *opt);

/**
 * Assemble a list of sequences
 *
 * @param opt      parameters
 * @param n_seqs   number of input sequences
 * @param seqs     sequences to assemble; FREED on return
 * @param n_utg    (out) number of unitigs in return
 *
 * @return array of unitigs
 */
fml_utg_t *fml_assemble(const fml_opt_t *opt, int n_seqs, bseq1_t *seqs, int *n_utg);

/**
 * Free unitigs
 *
 * @param n_utg    number of unitigs
 * @param utg      array of unitigs
 */
void fml_utg_destroy(int n_utg, fml_utg_t *utg);

/************************************************
 * Mid-level functions called by fml_assemble() *
 ************************************************/

/**
 * Adjust parameters based on input sequences
 *
 * @param opt       parameters to update IN PLACE
 * @param n_seqs    number of sequences
 * @param seqs      array of sequences
 */
void fml_opt_adjust(fml_opt_t *opt, int n_seqs, const bseq1_t *seqs);

/**
 * Error correction
 *
 * @param opt       parameters
 * @param n         number of sequences
 * @param seq       array of sequences; corrected IN PLACE
 *
 * @return k-mer coverage
 */
float fml_correct(const fml_opt_t *opt, int n, bseq1_t *seq);
float fml_fltuniq(const fml_opt_t *opt, int n, bseq1_t *seq);

/**
 * Construct FMD-index
 *
 * @param opt       parameters
 * @param n         number of sequences
 * @param seq       array of sequences; FREED on return
 *
 * @return FMD-index
 */
struct rld_t *fml_seq2fmi(const fml_opt_t *opt, int n, bseq1_t *seq);

/**
 * Generate initial overlap graph
 *
 * @param opt       parameters
 * @param e         FMD-index; FREED on return
 *
 * @return overlap graph in the "mag" structure
 */
struct mag_t *fml_fmi2mag(const fml_opt_t *opt, struct rld_t *e);

/**
 * Clean a mag graph
 *
 * @param opt       parameters
 * @param g         overlap graph; modified IN PLACE
 */
void fml_mag_clean(const fml_opt_t *opt, struct mag_t *g);

/**
 * Convert a graph in mag to fml_utg_t
 *
 * @param g         graph in the "mag" structure; FREED on return
 * @param n_utg     (out) number of unitigs
 *
 * @return array of unitigs
 */
fml_utg_t *fml_mag2utg(struct mag_t *g, int *n_utg);

/**
 * Output unitig graph in the mag format
 *
 * @param n_utg     number of unitigs
 * @param utg       array of unitigs
 */
void fml_utg_print(int n_utgs, const fml_utg_t *utg);

/**
 * Output unitig graph in the GFA format
 *
 * @param n_utg     number of unitigs
 * @param utg       array of unitigs
 */
void fml_utg_print_gfa(int n, const fml_utg_t *utg);

/**
 * Deallocate an FM-index
 *
 * @param e         pointer to the FM-index
 */
void fml_fmi_destroy(struct rld_t *e);

/**
 * Deallocate a mag graph
 *
 * @param g         pointer to the mag graph
 */
void fml_mag_destroy(struct mag_t *g);

#ifdef __cplusplus
}
#endif

#endif
