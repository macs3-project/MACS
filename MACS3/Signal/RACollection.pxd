cdef extern from "fml.h":
    ctypedef struct bseq1_t:
        int l_seq
        char *seq
        char *qual # NULL-terminated strings; length expected to match $l_seq

    ctypedef struct magopt_t:
        int flag, min_ovlp, min_elen, min_ensr, min_insr, max_bdist, max_bdiff, max_bvtx, min_merge_len, trim_len, trim_depth
        float min_dratio1, max_bcov, max_bfrac

    ctypedef struct fml_opt_t:
        int n_threads        # number of threads; don't use multi-threading for small data sets
        int ec_k             # k-mer length for error correction; 0 for auto estimate
        int min_cnt, max_cnt # both occ threshold in ec and tip threshold in cleaning lie in [min_cnt,max_cnt]
        int min_asm_ovlp     # min overlap length during assembly
        int min_merge_len    # during assembly, don't explicitly merge an overlap if shorter than this value
        magopt_t mag_opt     # graph cleaning options

    ctypedef struct fml_ovlp_t:
        unsigned int len_, from_, id_, to_ 
        #unit32_t from  # $from and $to: 0 meaning overlapping 5'-end; 1 overlapping 3'-end
        #unsigned int id
        #unsigned int to    # $id: unitig number

    ctypedef struct fml_utg_t:
        int len      # length of sequence
        int nsr      # number of supporting reads
        char *seq        # unitig sequence
        char *cov        # cov[i]-33 gives per-base coverage at i
        int n_ovlp[2]    # number of 5'-end [0] and 3'-end [1] overlaps
        fml_ovlp_t *ovlp # overlaps, of size n_ovlp[0]+n_ovlp[1]

    void fml_opt_init(fml_opt_t *opt)
    fml_utg_t* fml_assemble(const fml_opt_t *opt, int n_seqs, bseq1_t *seqs, int *n_utg)
    void fml_utg_destroy(int n_utg, fml_utg_t *utg)
    void fml_utg_print(int n_utgs, const fml_utg_t *utg)
    bseq1_t *bseq_read(const char *fn, int *n)

# --- end of fermi-lite functions ---

# --- smith-waterman alignment functions ---

cdef extern from "swalign.h":
    ctypedef struct seq_pair_t:
        char *a
        unsigned int alen
        char *b
        unsigned int blen
    ctypedef struct align_t:
        seq_pair_t *seqs
        char *markup;
        int start_a
        int start_b
        int end_a
        int end_b
        int matches
        int gaps
        double score
    align_t *smith_waterman(seq_pair_t *problem)
    void destroy_seq_pair(seq_pair_t *pair)
    void destroy_align(align_t *ali)
