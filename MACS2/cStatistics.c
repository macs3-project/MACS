//Time-stamp: <2014-07-29 17:01:45 Tao Liu>

#include <math.h>

float log10_poisson_cdf ( unsigned int n, float lam, short lower );
float log10_poisson_cdf_P_large_lambda ( unsigned int k, float lbd );
float log10_poisson_cdf_Q_large_lambda ( unsigned int k, float lbd );
float chi2_k1_cdf ( float x );
float chi2_k2_cdf ( float x );
float chi2_k4_cdf ( float x );
float log10_chi2_k1_cdf ( float x );
float log10_chi2_k2_cdf ( float x );
float log10_chi2_k4_cdf ( float x );

#define max(x, y) ((x)>(y)?(x):(y))
#define logspace_add(x, y) ((x)>(y)?(x):(y)) + log1p ( exp( -fabs(x - y) ) )
#define LOG10_E 0.43429448190325176

float chi2_k1_cdf ( float x )
{
  /*CDF for chi-square distribution with degree of freedom 1.*/
  return erf( sqrt(x/2) );
}

float log10_chi2_k1_cdf ( float x )
{
  /*log10 CDF for chi-square distribution with degree of freedom 1.*/
  return log10( erf( sqrt(x/2) ) );
}

float chi2_k2_cdf ( float x )
{
  /*CDF for chi-square distribution with degree of freedom 2.*/
  return 1 - exp( -x/2 );
}

float log10_chi2_k2_cdf ( float x )
{
  /*log10 CDF for chi-square distribution with degree of freedom 2.*/
  return log1p( - exp( -x/2 ) ) * LOG10_E;
}

float chi2_k4_cdf ( float x )
{
  /*CDF for chi-square distribution with degree of freedom 4.*/
  return 1 - exp( -x/2 ) * ( 1 + x/2 );
}

float log10_chi2_k4_CDF ( float x )
{
  /*log10 CDF for chi-square distribution with degree of freedom 4.*/
  return log1p( - exp( -x/2 ) * ( 1 + x/2 ) ) * LOG10_E;
}


float log10_poisson_cdf ( unsigned int n, float lam, short lower )
{
  /*Poisson CDF evaluater.
    
    This is a more stable CDF function. It can tolerate large lambda
    value. While the lambda is larger than 700, the function will be a
    little slower.

    Parameters:
    n     : your observation
    lam   : lambda of poisson distribution
    lower : if lower is False, calculate the upper tail CDF, otherwise, to calculate lower tail; Default is False.
    log10 : if log10 is True, calculation will be in log space. Default is False.
  */
  if ( lower )
    // lower tail
    return log10_poisson_cdf_P_large_lambda( n, lam );
  else
    // upper tail
    return log10_poisson_cdf_Q_large_lambda( n, lam );
}

float log10_poisson_cdf_P_large_lambda ( unsigned int k, float lbd )
{
  /*Slower Poisson CDF evaluater for lower tail which allow
    calculation in log space. Better for the pvalue < 10^-310.

    Parameters:
    k	: observation
    lbd	: lambda

    ret = -lambda + \ln( \sum_{i=k+1}^{\inf} {lambda^i/i!} = -lambda + \ln( sum{ exp{ln(F)} } ), where F=lambda^m/m!
    \ln{F(m)} = m*ln{lambda} - \sum_{x=1}^{m}\ln(x)
    Calculate \ln( sum{exp{N} ) by logspace_add function

    Return the log10(pvalue)
  */
  float residue = 0;
  float logx = 0;
  float ln_lbd = log( lbd );
  float logy, pre_residue;

  // first residue
  int m = k;
  float sum_ln_m = 0;
  int i = 0;
  for ( i = 1; i <= m; i++ )
    sum_ln_m += log( i );
  logx = m * ln_lbd - sum_ln_m;
  residue = logx;

  for ( ; m > 1; m-- )
    {
      logy = logx - ln_lbd + log( m );
      pre_residue = residue;
      residue = logspace_add( pre_residue, logy );
      if ( fabs( pre_residue - residue ) < 1e-5 ) break;
      logx = logy;
    }
  
  return ( residue-lbd ) / M_LN10 ;
}

float log10_poisson_cdf_Q_large_lambda ( unsigned int k, float lbd )
{
  /*Slower Poisson CDF evaluater for upper tail which allow
    calculation in log space. Better for the pvalue < 10^-310.

    Parameters:
    k	: observation
    lbd	: lambda

    ret = -lambda + \ln( \sum_{i=k+1}^{\inf} {lambda^i/i!} = -lambda + \ln( sum{ exp{ln(F)} } ), where F=lambda^m/m!
    \ln{F(m)} = m*ln{lambda} - \sum_{x=1}^{m}\ln(x)
    Calculate \ln( sum{exp{N}} ) by logspace_add function

    Return the log10(pvalue)
  */
  float residue = 0;		/* to save ln(\sum{ exp{N} }) */
  float logx = 0;
  float ln_lbd = log( lbd );
  float logy, pre_residue;

  // first residue
  int m = k + 1 ;
  float sum_ln_m = 0;
  int i;

  for ( i = 2; i <= m; i++ )	/* should start from 1, however log1 = 0 */
    sum_ln_m += log( i );
  logx = m * ln_lbd - sum_ln_m;
  residue = logx;		/* first residue calculated here. */

  while ( 1 )
    {
      m++;
      logy = logx + ln_lbd - log( m );
      pre_residue = residue;
      residue = logspace_add( pre_residue, logy );
      if ( fabs( pre_residue - residue ) < 1e-5 ) break;
      logx = logy;
    }
  
  return (residue - lbd ) / M_LN10; 
}
