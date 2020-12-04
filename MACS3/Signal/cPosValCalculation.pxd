# cython: language_level=3

# ------------------------------------
# NumPy modules
# ------------------------------------
import numpy as np
cimport numpy as np
from numpy cimport int32_t
ctypedef np.float32_t float32_t

cdef extern from "cPosValCalculation.h":
    cdef struct PosVal:
        int pos
        float value

    PosVal * single_end_pileup ( int * plus_tags, long l_plus_tags, int * minus_tags, long l_minus_tags, int five_shift, int three_shift, int leftmost_coord, int rightmost_coord,  float scale_factor, float baseline_value, long * final_length ) nogil

    PosVal * quick_pileup ( int * start_poss, int * end_poss, long length_poss, float scale_factor, float baseline_value, long * final_length ) nogil

    int * fix_coordinates ( int * poss, long l, int leftmost_coord, int rightmost_coord ) nogil

    PosVal * max_over_two_pv_array ( PosVal * pva1, long l_pva1, PosVal * pva2, long l_pva2, long * final_length ) nogil

    void write_pv_array_to_bedGraph ( PosVal * pv_array, long l_pv_array, char * chromosome, char * bdgfile, short append )

    long quick_pileup_simple ( int32_t * ret_poss, float32_t * ret_values, int32_t * start_poss, int32_t * end_poss, long length_poss, float scale_factor, float baseline_value )

