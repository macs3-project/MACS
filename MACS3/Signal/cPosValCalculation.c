/* Time-stamp: <2020-11-23 14:52:31 Tao Liu>

 This code is free software; you can redistribute it and/or modify it
 under the terms of the BSD License (see the file LICENSE included
 with the distribution).
*/

#include <stdio.h>
#include <stdlib.h>

#include "cPosValCalculation.h"

/* Simple compare function */
int cmpfunc_simple (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

/* Fix coordinates

   Input:

   1. int * poss: must be sorted!

   Return:

*/
int * fix_coordinates ( int * poss, long l, int leftmost_coord, int rightmost_coord )
{
  long i;
  for ( i = 0; i < l; i++ )
    {
      if ( poss[ i ] < leftmost_coord )
	poss[ i ] = leftmost_coord;
      else
	break;
    }

  for ( i = l-1; i > -1; i-- )
    {
      if ( poss[ i ] > rightmost_coord )
	poss[ i ] = rightmost_coord;
      else
	break;
    }
  return poss;
}

/* Pileup function for single end data.
   Input:
   1. int * plus_tags: the 5' ends of tags aligned to plus strand. Start from 0.
   2. long l_plus_tags:  number of plus tags.
   3. int * minus_tags: the 3' ends of tags aligned to minus strand. Start from 0.
   4. long l_minus_tags: number of minus tags.
   5. int five_shift: shift value at 5' end to recover the DNA fragment
   6. int three_shift: shift value at 3' end to recover the DNA fragment
   7. int leftmost_coord: any coords smaller than it will be set to it
   8. int rightmost_coord: any coords larger than it will be set to it
   9. float scale_factor: scale factor on pileup
   10. float baseline_value: the minimum value of pileup

   Return:
   1. PosVal *: position-value pairs in bedGraph fashion.
   2. long * final_length: the final length for the returned array of PosVal.

   Example of usage:

   pileup = single_end_pileup ( {1,2,3}, 3, {2,3,4}, 3, 0, 3, 0, 20, 0.5, 1.0, &final_length );
   for ( i = 0; i < final_length; i++ )
    {
      printf( "pos:%d value:%.2f\n", pileup[i].pos, pileup[i].value );
    }

 */
struct PosVal * single_end_pileup( int * plus_tags, long l_plus_tags, int * minus_tags, long l_minus_tags, int five_shift, int three_shift, int leftmost_coord, int rightmost_coord, float scale_factor, float baseline_value, long * final_length )
{
  long i_s, i_e, i;
  int p, pre_p, pileup;
  long l = l_plus_tags + l_minus_tags;
  int * start_poss, * end_poss, * ptr_start_poss, * ptr_end_poss;
  struct PosVal * pos_value_array;

  start_poss = (int *) malloc( l * sizeof( int ) );
  end_poss = (int *) malloc( l * sizeof( int ) );

  ptr_start_poss = start_poss;
  ptr_end_poss = end_poss;

  for ( i = 0; i < l_plus_tags; i++ )
    {
      *ptr_start_poss = plus_tags[ i ] - five_shift; ptr_start_poss++;
      *ptr_end_poss   = plus_tags[ i ] + three_shift; ptr_end_poss++;
    }

  for ( i = 0; i < l_minus_tags; i++ )
    {
      *ptr_start_poss = minus_tags[ i ] - three_shift; ptr_start_poss++;
      *ptr_end_poss   = minus_tags[ i ] + five_shift; ptr_end_poss++;
    }

  qsort( start_poss, l, sizeof( int ), cmpfunc_simple );
  qsort( end_poss,   l, sizeof( int ), cmpfunc_simple );

  // fix negative coordinations and those extends over end of chromosomes
  start_poss = fix_coordinates ( start_poss, l,  leftmost_coord, rightmost_coord );
  end_poss   = fix_coordinates ( end_poss, l, leftmost_coord, rightmost_coord );

  pos_value_array = quick_pileup ( start_poss, end_poss, l, scale_factor, baseline_value, final_length );

  // clean mem
  free( start_poss );
  free( end_poss );
  return pos_value_array;
}

/* Assume start_poss and end_poss have been sorted and the coordinates have been fixed. */
struct PosVal * quick_pileup ( int * start_poss, int * end_poss, long length_poss, float scale_factor, float baseline_value, long * final_length )
{
  long i_s, i_e, i, I;
  int p, pre_p, pileup;
  struct PosVal * pos_value_array, * ptr_pos_value_array;
  int * ptr_start_poss, * ptr_end_poss;
  long l = length_poss;

  pos_value_array = (struct PosVal *) malloc ( 2 * l * sizeof( struct PosVal ) );

  i_s = 0; i_e = 0;

  ptr_pos_value_array = pos_value_array; ptr_start_poss = start_poss; ptr_end_poss = end_poss;

  pileup = 0;
  pre_p = min( *start_poss, *end_poss );
  ptr_start_poss++; ptr_end_poss++;

  I = 0;
  if ( pre_p != 0 )
    {
      (*ptr_pos_value_array).pos = pre_p;
      (*ptr_pos_value_array).value = max( 0, baseline_value );
      ptr_pos_value_array++; I++;
    }

  ptr_start_poss = start_poss;
  ptr_end_poss = end_poss;

  while (i_s < l && i_e < l)
    {
      if ( *ptr_start_poss < *ptr_end_poss )
	{
	  p = *ptr_start_poss;
	  if ( p != pre_p )
	    {
	      (*ptr_pos_value_array).pos = p;
	      (*ptr_pos_value_array).value = max( pileup * scale_factor, baseline_value );
	      ptr_pos_value_array++; I++;
	      pre_p = p;
	    }
	  pileup += 1;
	  i_s += 1;
	  ptr_start_poss++;
	}
      else if ( *ptr_start_poss > *ptr_end_poss )
	{
	  p = *ptr_end_poss;
	  if ( p != pre_p )
	    {
	      (*ptr_pos_value_array).pos = p;
	      (*ptr_pos_value_array).value = max( pileup * scale_factor, baseline_value );
	      ptr_pos_value_array++; I++;
	      pre_p = p;
	    }
	  pileup -= 1;
	  i_e += 1;
	  ptr_end_poss++;
	}
      else
	{
	  i_s += 1;
	  i_e += 1;
	  ptr_start_poss++;
          ptr_end_poss++;
	}
    }

  // add the rest of end positions.
  if ( i_e < l )
    {
      for ( i = i_e; i < l; i++ )
	{
	  p = *ptr_end_poss;
	  if ( p != pre_p )
	    {
	      (*ptr_pos_value_array).pos = p;
	      (*ptr_pos_value_array).value = max( pileup * scale_factor, baseline_value );
	      ptr_pos_value_array++; I++;
	      pre_p = p;
	    }
	  pileup -= 1;
	  ptr_end_poss++;
	}
    }
  pos_value_array = (struct PosVal *) realloc ( pos_value_array, I * sizeof( struct PosVal ) );
  *final_length = I;		/* return the final length of pos_value_array */
  return pos_value_array;
}

long quick_pileup_simple ( int * ret_poss, float * ret_values, int * start_poss, int * end_poss, long length_poss, float scale_factor, float baseline_value )
{
  long i_s, i_e, i, I;
  int p, pre_p, pileup;
  int * ptr_ret_poss;
  int * ptr_start_poss, * ptr_end_poss;
  float * ptr_ret_values;
  long l = length_poss;

  ptr_ret_poss = ret_poss; ptr_ret_values = ret_values;
  ptr_start_poss = start_poss; ptr_end_poss = end_poss;

  i_s = 0; i_e = 0;

  pileup = 0;
  pre_p = min( *ptr_start_poss, *ptr_end_poss );
  ptr_start_poss++; ptr_end_poss++;

  I = 0;
  if ( pre_p != 0 )
    {
      *ptr_ret_poss = pre_p;
      *ptr_ret_values = max( 0, baseline_value );
      ptr_ret_poss++; ptr_ret_values++; I++;
    }

  while (i_s < l && i_e < l)
    {
      if ( *ptr_start_poss < *ptr_end_poss )
	{
	  p = *ptr_start_poss;
	  if ( p != pre_p )
	    {
	      *ptr_ret_poss = p;
	      *ptr_ret_values = max( pileup * scale_factor, baseline_value );
	      ptr_ret_poss++;ptr_ret_values++; I++;
	      pre_p = p;
	    }
	  pileup += 1;
	  i_s += 1;
	  ptr_start_poss++;
	}
      else if ( *ptr_start_poss > *ptr_end_poss )
	{
	  p = *ptr_end_poss;
	  if ( p != pre_p )
	    {
	      *ptr_ret_poss = p;
	      *ptr_ret_values = max( pileup * scale_factor, baseline_value );
	      ptr_ret_poss++;ptr_ret_values++; I++;
	      pre_p = p;
	    }
	  pileup -= 1;
	  i_e += 1;
	  ptr_end_poss++;
	}
      else
	{
	  i_s += 1;
	  i_e += 1;
	  ptr_start_poss++;
          ptr_end_poss++;
	}
    }

  // add the rest of end positions.
  if ( i_e < l )
    {
      for ( i = i_e; i < l; i++ )
	{
	  p = *ptr_end_poss;
	  if ( p != pre_p )
	    {
	      *ptr_ret_poss = p;
	      *ptr_ret_values = max( pileup * scale_factor, baseline_value );
	      ptr_ret_poss++;ptr_ret_values++; I++;
	      pre_p = p;
	    }
	  pileup -= 1;
	  ptr_end_poss++;
	}
    }
  return I;
}

/* Calculate the maximum value between two sets of PosVal arrays (like bedGraph type of data) */
struct PosVal * max_over_two_pv_array ( struct PosVal * pva1, long l_pva1, struct PosVal * pva2, long l_pva2, long * final_length )
{
  struct PosVal * ptr_pva1, * ptr_pva2;
  struct PosVal * ret_pva, * ptr_ret_pva;
  long i, i1, i2, I;

  ptr_pva1 = pva1; ptr_pva2 = pva2;
  ret_pva = ( struct PosVal * ) malloc ( ( l_pva1 + l_pva2 ) * sizeof( struct PosVal ) );
  ptr_ret_pva = ret_pva;

  i1 = i2 = 0;

  I = 0;

  while ( i1< l_pva1 && i2 < l_pva2 )
    {
      if ( (*ptr_pva1).pos < (*ptr_pva2).pos )
	{
	  (*ptr_ret_pva).pos = (*ptr_pva1).pos;
	  (*ptr_ret_pva).value = max( (*ptr_pva1).value, (*ptr_pva2).value );
	  ptr_ret_pva++;I++;
	  ptr_pva1++; i1++;
	}
      else if ( (*ptr_pva1).pos > (*ptr_pva2).pos )
	{
	  (*ptr_ret_pva).pos = (*ptr_pva2).pos;
	  (*ptr_ret_pva).value = max( (*ptr_pva1).value, (*ptr_pva2).value );
	  ptr_ret_pva++;I++;
	  ptr_pva2++; i2++;
	}
      else // (*ptr_pva1).pos == (*ptr_pva2).pos
	{
	  (*ptr_ret_pva).pos = (*ptr_pva1).pos;
	  (*ptr_ret_pva).value = max( (*ptr_pva1).value, (*ptr_pva2).value );
	  ptr_ret_pva++;I++;
	  ptr_pva1++; i1++;
	  ptr_pva2++; i2++;
	}
    }

  *final_length = I;
  return ret_pva;
}

/* Calculate using specified function between two sets of PosVal arrays (like bedGraph type of data) */
struct PosVal * apply_func_two_pv_array ( float (*func)(float, float), struct PosVal * pva1, long l_pva1, struct PosVal * pva2, long l_pva2, long * final_length )
{
  struct PosVal * ptr_pva1, * ptr_pva2;
  struct PosVal * ret_pva, * ptr_ret_pva;
  long i, i1, i2, I;

  ptr_pva1 = pva1; ptr_pva2 = pva2;
  ret_pva = ( struct PosVal * ) malloc ( ( l_pva1 + l_pva2 ) * sizeof( struct PosVal ) );
  ptr_ret_pva = ret_pva;

  i1 = i2 = 0;

  I = 0;

  while ( i1< l_pva1 && i2 < l_pva2 )
    {
      if ( (*ptr_pva1).pos < (*ptr_pva2).pos )
	{
	  (*ptr_ret_pva).pos = (*ptr_pva1).pos;
	  (*ptr_ret_pva).value = (*func)( (*ptr_pva1).value, (*ptr_pva2).value );
	  ptr_ret_pva++;I++;
	  ptr_pva1++; i1++;
	}
      else if ( (*ptr_pva1).pos > (*ptr_pva2).pos )
	{
	  (*ptr_ret_pva).pos = (*ptr_pva2).pos;
	  (*ptr_ret_pva).value = (*func)( (*ptr_pva1).value, (*ptr_pva2).value );
	  ptr_ret_pva++;I++;
	  ptr_pva2++; i2++;
	}
      else // (*ptr_pva1).pos == (*ptr_pva2).pos
	{
	  (*ptr_ret_pva).pos = (*ptr_pva1).pos;
	  (*ptr_ret_pva).value = (*func)( (*ptr_pva1).value, (*ptr_pva2).value );
	  ptr_ret_pva++;I++;
	  ptr_pva1++; i1++;
	  ptr_pva2++; i2++;
	}
    }

  *final_length = I;
  return ret_pva;
}

/* Align two PosVal arrays according to their overlaps */
struct PosValVal * align_two_pv_array ( struct PosVal * pva1, long l_pva1, struct PosVal * pva2, long l_pva2, long * final_length )
{
  struct PosVal * ptr_pva1, * ptr_pva2;
  struct PosValVal * ret_pvva, * ptr_ret_pvva;
  long i, i1, i2, I;

  ptr_pva1 = pva1; ptr_pva2 = pva2;
  ret_pvva = ( struct PosValVal * ) malloc ( ( l_pva1 + l_pva2 ) * sizeof( struct PosValVal ) );
  ptr_ret_pvva = ret_pvva;

  i1 = i2 = 0;

  I = 0;

  while ( i1< l_pva1 && i2 < l_pva2 )
    {
      if ( (*ptr_pva1).pos < (*ptr_pva2).pos )
	{
	  (*ptr_ret_pvva).pos = (*ptr_pva1).pos;
	  (*ptr_ret_pvva).value1 = (*ptr_pva1).value;
	  (*ptr_ret_pvva).value2 = (*ptr_pva2).value;
	  ptr_ret_pvva++;I++;
	  ptr_pva1++; i1++;
	}
      else if ( (*ptr_pva1).pos > (*ptr_pva2).pos )
	{
	  (*ptr_ret_pvva).pos = (*ptr_pva2).pos;
	  (*ptr_ret_pvva).value1 = (*ptr_pva1).value;
	  (*ptr_ret_pvva).value2 = (*ptr_pva2).value;
	  ptr_ret_pvva++;I++;
	  ptr_pva2++; i2++;
	}
      else // (*ptr_pva1).pos == (*ptr_pva2).pos
	{
	  (*ptr_ret_pvva).pos = (*ptr_pva1).pos;
	  (*ptr_ret_pvva).value1 = (*ptr_pva1).value;
	  (*ptr_ret_pvva).value2 = (*ptr_pva2).value;
	  ptr_ret_pvva++;I++;
	  ptr_pva1++; i1++;
	  ptr_pva2++; i2++;
	}
    }

  *final_length = I;
  return ret_pvva;
}

/* Write pos-value array to a bedGraph file. If append is non-zero then just add content to the existing file. */
void write_pv_array_to_bedGraph ( struct PosVal * pv_array, long l_pv_array, char * chromosome, char * bdgfile, short append )
{
  int pre_s, pre_e;
  float pre_v;
  long i;
  FILE * fp;

  if ( append > 0 )
    fp = fopen ( bdgfile, "a" );
  else
    fp = fopen ( bdgfile, "w" );

  pre_s = 0;

  pre_e = (*pv_array).pos;
  pre_v = (*pv_array).value;

  pv_array += 1;

  for ( i = 1; i < l_pv_array; i++ )
    {
      if ( (*pv_array).value != pre_v )
	{
	  fprintf ( fp, "%s\t%d\t%d\t%.5f\n", chromosome, pre_s, pre_e, pre_v );
	  pre_s = pre_e;
	  pre_e = (*pv_array).pos;
	  pre_v = (*pv_array).value;
	}
      else
	{
	  pre_e = (*pv_array).pos;
	}
      pv_array ++;
    }

  /* last piece */
  fprintf ( fp, "%s\t%d\t%d\t%.5f\n", chromosome, pre_s, pre_e, pre_v );

  fclose( fp );
}

/* for testing */
int main()
{
  int i;
  int five_shift = 0;
  int three_shift = 2;
  int p_array1[3] = {11,12,13};
  int p_array2[3] = {12,13,14};
  int m_array1[3] = {12,13,14};
  int m_array2[3] = {13,14,15};
  int leftmost_coord = 0;
  int rightmost_coord = 20;
  float scale_factor = 0.5;
  long final_length1 = 0;
  long final_length2 = 0;
  long final_length_max = 0;
  struct PosVal * pileup1;
  struct PosVal * pileup2;
  struct PosVal * max_pileup;

  pileup1 = single_end_pileup ( p_array1, 3, m_array1, 3, five_shift, three_shift, leftmost_coord, rightmost_coord, scale_factor, 0, &final_length1 );
  pileup2 = single_end_pileup ( p_array2, 3, m_array2, 3, five_shift, three_shift, leftmost_coord, rightmost_coord, scale_factor, 0, &final_length2 );

  printf( "pileup 1\n" );
  for ( i = 0; i < final_length1; i++ )
    {
      printf( "pos:%d value:%.2f\n", pileup1[i].pos, pileup1[i].value );
    }
  printf( "pileup 2\n" );
  for ( i = 0; i < final_length2; i++ )
    {
      printf( "pos:%d value:%.2f\n", pileup2[i].pos, pileup2[i].value );
    }

  max_pileup = max_over_two_pv_array ( pileup1, final_length1, pileup2, final_length2, &final_length_max );
  printf( "max of pileup 1 and 2\n" );
  for ( i = 0; i < final_length_max; i++ )
    {
      printf( "pos:%d value:%.2f\n", max_pileup[i].pos, max_pileup[i].value );
    }

  printf( "write to bedGraph\n" );

  write_pv_array_to_bedGraph ( max_pileup, final_length_max, "chr1", "test.bdg", 0 );

}



