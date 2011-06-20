from cpython cimport bool

cpdef factorial ( unsigned int n )

cpdef poisson_cdf ( unsigned int n, double lam, bool lower = ?, bool log10 = ? )

cdef __poisson_cdf ( unsigned int k, double a )
    
cdef __poisson_cdf_large_lambda ( unsigned int k, double a )

cdef __poisson_cdf_Q ( unsigned int k, double a )

cdef __poisson_cdf_Q_large_lambda ( unsigned int k, double a )

cdef double log10_poisson_cdf_P_large_lambda ( unsigned int k, double lbd )

cdef double log10_poisson_cdf_Q_large_lambda ( unsigned int k, double lbd )

cdef double logspace_add ( double logx, double logy )

cpdef poisson_cdf_inv ( double cdf, double lam, int maximum = ? )

cpdef poisson_cdf_Q_inv ( double cdf, double lam, int maximum = ? )

cpdef poisson_pdf ( unsigned int k, double a )

cdef binomial_coef ( long n, long k )

cpdef binomial_cdf ( long x, long a, double b, bool lower = ? )

cdef _binomial_cdf_r ( long x, long a, double b )

cdef _binomial_cdf_f ( long x, long a, double b )

cpdef binomial_cdf_inv ( double cdf, long a, double b )
    
cpdef binomial_pdf( long x, long a, double b )
