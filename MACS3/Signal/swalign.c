/*  
 *  Copyright (c) 2010 Nicolaus Lance Hepler
 *  
 * 
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use,
 *  copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following
 *  conditions:
 * 
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 * 
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *  OTHER DEALINGS IN THE SOFTWARE.
 */
// Note: That's an MIT license.
// All double-commented comments below are from Nicolaus Lance Hepler.
// Original repository: https://code.google.com/archive/p/swalign/

/* Tao Liu made some modifications.... Note to myself: maybe I should use ksw.c/h instead */

#include "swalign.h"

#define GAPO -10.0
#define GAPE -2.0
#define MATCH 2.0
#define MISMATCH -3.0		/* the scoring parameters are from BLASTN default */
//             ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz
#define TRANS "TVGHEFCDIJMLKNOPQYWAABSXRZ[\\]^_`tvghefcdijmlknopqywaabsxrz"
#define TRANS_OFFSET 65

#define STOP 0
#define LEFT 1
#define DIAGONAL 2
#define UP 3

// /* reverse a string str0 in place, return str */
static char* reverse(char *str, unsigned int l) {
  char *left  = str;
  char *right = left + l - 1;
  char tmp;
  while (left < right) {
    tmp        = *left;
    *(left++)  = *right;
    *(right--) = tmp;
  }
  *(str+l) = '\0';
  return str;
}

// Return the reverse complement of a sequence.
char* revcomp(char *str) {
  char *left = str;
  char *right = left + strlen(str) - 1;
  char tmp;

  while (left < right) {
    tmp        = get_char_comp(*left);
    *(left++)  = get_char_comp(*right);
    *(right--) = tmp;
  }

  return str;
}

// Return the complement of a base.
// Uses a simple lookup table: a string with the complements of all possible sequence characters.
static char get_char_comp(char c) {
  int i = c - TRANS_OFFSET;
  if (i < 0 || i > 57) {
    return c;
  } else {
    return TRANS[i];
  }
}

// return the alignment of two sequences
static align_t *traceback(seq_pair_t *problem, unsigned short *S, int row, int col, float score, unsigned short *ngap_vertical, unsigned short *ngap_horizontal) {
  unsigned int l, l2;
  unsigned int m = problem->alen + 1;
  unsigned int n = problem->blen + 1;
  align_t *result = malloc(sizeof(align_t));
  seq_pair_t *seqs = malloc(sizeof(seq_pair_t));

  int max_len = problem->alen + problem->blen; /* maximum length of alignment */

  char reversed1[max_len];	/* reversed sequence #1 */
  char reversed2[max_len];	/* reversed sequence #2 */
  char reversed3[max_len];	/* reversed markup */

  unsigned int len1 = 0;	/* length of seq #1 in alignment */
  unsigned int len2 = 0;	/* length of seq #2 in alignment */
  unsigned int len3 = 0;	/* length of the markup line in alignment */  

  unsigned int identity = 0; // count of identitcal pairs
  unsigned int gaps = 0; // count of gaps

  char c1, c2;

  int i = row; // traceback start row
  int j = col; // traceback start col
  int k = i * n;
  bool still_going = true; // traceback flag: true -> continue & false -> stop

  while ( still_going ) {
    switch ( S[k+j] ) {
    case UP:
      for (l = 0, l2 = ngap_vertical[k + j]; l < l2; l++) {
	reversed1[len1++] = problem->a[--i];
	reversed2[len2++] = '-';
	reversed3[len3++] = ' ';
	k -= n;
	gaps++;
      }
      break;
    case DIAGONAL:
      c1 = problem->a[--i];
      c2 = problem->b[--j];
      k -= n;
      reversed1[len1++] = c1;
      reversed2[len2++] = c2;
      if (c1 == c2) {
	reversed3[len3++] = '|';
	identity++;
      } else
	reversed3[len3++] = '.';
      break;
    case LEFT:
      for (l = 0, l2 = ngap_horizontal[k + j]; l < l2; l++) {
	reversed1[len1++] = '-';
	reversed2[len2++] = problem->b[--j];
	reversed3[len3++] = ' ';
	gaps++;
      }
      break;
    case STOP:
      still_going = false;
    }
  }

  seqs->a = malloc(sizeof(char) * (len1 + 1) );
  seqs->b = malloc(sizeof(char) * (len2 + 1) );
  result->markup = malloc(sizeof(char) * (len3 + 1) ); 
  
  memset(seqs->a, '\0', sizeof(char) * (len1 + 1));
  memset(seqs->b, '\0', sizeof(char) * (len2 + 1));
  memset(result->markup, '\0', sizeof(char) * (len3 + 1) );

  reverse(reversed1, len1);
  reverse(reversed2, len2);
  reverse(reversed3, len3);

  strcpy(seqs->a, reversed1);
  strcpy(seqs->b, reversed2);
  strcpy(result->markup, reversed3);

  seqs->alen = k;
  seqs->blen = k;

  result->seqs = seqs;
  result->score = score;
  result->matches = identity;
  result->gaps = gaps;
  result->start_a = i;
  result->start_b = j;
  result->end_a = row;
  result->end_b = col;

  //printf("%d %d %d",len1, len2, len3);
  
  return result; 
}

void destroy_seq_pair(seq_pair_t *pair) {
  free(pair->a);
  free(pair->b);
  free(pair);
  return;
}

void destroy_align(align_t *ali) {
  destroy_seq_pair( ali->seqs );
  free(ali->markup);
  return;
}

//align_t *smith_waterman(seq_pair_t *problem, bool local) {
align_t *smith_waterman(seq_pair_t *problem) {
  unsigned int m = problem->alen + 1;
  unsigned int n = problem->blen + 1;

  /* traceback matrix */
  unsigned short *S = malloc(sizeof(unsigned short) * m * n); /* 0 = STOP; 1 = LEFT; 2 = DIAGONAL; 3 = UP */

  /* number of vertical gaps for each cell */
  unsigned short *ngap_vertical = malloc(sizeof(unsigned short) * m * n);
  /* number of horizontal gaps for each cell */
  unsigned short *ngap_horizontal = malloc(sizeof(unsigned short) * m * n);

  align_t *result;
  unsigned int i, j, k, l;

  float f;			         /* score of final alignment */
  float * g = malloc(sizeof(float) * n); /* score if x_i aligns to a gap after y_i */
  float h;			         /* score if y_i aligns to a gap after x_i */
  float * v = malloc(sizeof(float) * n); /* best of score of alignment x_1 ... x_i to y_1 ... y_i */
  float v_diagonal;

  float sim_score, g1, g2, h1, h2;

  /* row and col number, and score of optimal cell in the matrix */
  unsigned int t_row, t_col;
  float t_score;

  /* Initialize traceback matrix */
  for ( i = 0, k = 0; i < m; i++, k+= n )
    S[k] = STOP;
  for ( j = 1; j < n; j++ )
    S[j] = STOP;

  /* set number of gaps */
  for (i = 0, k = 0; i < m; i++, k += n)
    for (j = 0; j < n; j++)
      ngap_vertical[k + j] = ngap_horizontal[k + j] = 1;
  
  
  g[0] = -INFINITY;
  h = -INFINITY;
  v[0] = 0;

  t_row = 0;
  t_col = 0;
  t_score = -INFINITY;

  for ( j = 1; j < n; j++ ) {
    g[j] = -INFINITY;
    v[j] = 0;
  }
  
  for (i = 1, k = n; i < m; i++, k += n) { /* i for row number/# on seq1; k for position of first column in m*n matrix; */
    h = -INFINITY;
    v_diagonal = v[0];
    for (j = 1, l = k + 1; j < n; j++, l++) { /* j for column number/# on seq2; l for position of ith row jth column in m*n matrix; */
      sim_score = (strncmp(problem->a+(i-1), problem->b+(j-1), 1) == 0) ? MATCH : MISMATCH;

      /* set direction in traceback matrix */
      f = v_diagonal + sim_score;

      g1 = g[j] + GAPE;
      g2 = v[j] + GAPO;

      if ( g1 > g2 ) {		/* perfer gap extension in vertical direction (seq 1) */
	g[j] = g1;
	ngap_vertical[l] = (short) (ngap_vertical[l - n] + 1);
      } else {			/* prefer gap openning */
	g[j] = g2;
      }

      h1 = h + GAPE;
      h2 = v[j-1] + GAPO;

      if (h1 > h2) {		/* prefer gap extention in horizontal direction ( seq 2) */
	h = h1;
	ngap_horizontal[l] = (short) (ngap_horizontal[l - 1] + 1);
      } else {
	h = h2;
      }
	  
      v_diagonal = v[j];
      //v[j] = max( f, g[j], h, 0 ); /* the final score */

      /* Determine the traceback direction */
      if ( f <= 0 && g[j] <=0 && h <=0 ) { /* 0 is max */
	  v[j] = 0;
	  S[l] = STOP;
	}
      else if ( g[j] <= f && h <= f ) { /* or f is max */
	v[j] = f;
	S[l] = DIAGONAL;
      }
      else if ( h <= g[j]) { 	/* or g[j] is max */
	v[j] = g[j];
	S[l] = UP;
      }
      else {			/* or h is max */
	v[j] = h;
	S[l] = LEFT;
      }
      
      // Set the traceback start at the current cell i, j and score
      if (v[j] > t_score) {
	t_row = i;
	t_col = j;
	t_score = v[j];
      }
    }
  }

  result = traceback(problem, S, t_row, t_col, t_score, ngap_vertical, ngap_horizontal);

  // print_matrix(S, problem);

  free(S);
  free(g);
  free(v);
  free(ngap_vertical);
  free(ngap_horizontal);
  
  return result;
}

void print_alignment(align_t *result) {
  printf("Score: %0.0f  Matches: %d Gaps: %d\n", result->score, result->matches, result->gaps);
  printf("Target: %3d %s %-3d\n", result->start_a, result->seqs->a, result->end_a);
  printf("            %s     \n", result->markup);
  printf("Query:  %3d %s %-3d\n", result->start_b, result->seqs->b, result->end_b);
}

int main(int argc, const char **argv) {

  if (argc != 3) {
    printf("usage: swalign TARGET_SEQ QUERY_SEQ\n");
    exit(1);
  }

  {
    seq_pair_t problem;
    align_t *result;
    char c[strlen(argv[1])], d[strlen(argv[2])];
  
    strcpy(c, argv[1]);
    strcpy(d, argv[2]);
  
    problem.a = c;
    problem.alen = strlen(problem.a);
    problem.b = d;
    problem.blen = strlen(problem.b);
  
    result = smith_waterman(&problem);
  
    print_alignment(result);
  }

  exit(0);
} 



