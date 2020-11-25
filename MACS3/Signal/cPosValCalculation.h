/* Time-stamp: <2020-11-20 16:37:31 Tao Liu>

 This code is free software; you can redistribute it and/or modify it
 under the terms of the BSD License (see the file LICENSE included
 with the distribution).
*/

#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)>(y)?(y):(x))
#define compare(x, y) ((x)>(y)?1:0)

/* Equivalent to bedGraph data structure */
struct PosVal {
  int pos;
  float value;
};

/* for comparison between two PosVal arrays */
struct PosValVal {
  int pos;
  float value1;
  float value2;
};

struct PosVal * single_end_pileup ( int * plus_tags, long l_plus_tags, int * minus_tags, long l_minus_tags, int five_shift, int three_shift, int leftmost_coord, int rightmost_coord,  float scale_factor, float baseline_value, long * final_length );

struct PosVal * quick_pileup ( int * start_poss, int * end_poss, long length_poss, float scale_factor, float baseline_value, long * final_length );

int cmpfunc_simple ( const void * a, const void * b);

int * fix_coordinates ( int * poss, long l, int leftmost_coord, int rightmost_coord );

struct PosVal * max_over_two_pv_array ( struct PosVal * pva1, long l_pva1, struct PosVal * pva2, long l_pva2, long * final_length );

struct PosVal * apply_func_two_pv_array ( float (*func)(float, float), struct PosVal * pva1, long l_pva1, struct PosVal * pva2, long l_pva2, long * final_length );

struct PosValVal * align_two_pv_array ( struct PosVal * pva1, long l_pva1, struct PosVal * pva2, long l_pva2, long * final_length );

void write_pv_array_to_bedGraph ( struct PosVal * pv_array, long l_pv_array, char * chromosome, char * bdgfile, short append );

long quick_pileup_simple ( int * ret_poss, float * ret_values, int * start_poss, int * end_poss, long length_poss, float scale_factor, float baseline_value );
