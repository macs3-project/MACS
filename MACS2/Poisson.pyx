cdef extern from "cPoisson.h":
    double log10_poisson_cdf ( unsigned int n, double lam, short lower )
    double log10_poisson_cdf_P_large_lambda ( unsigned int k, double lbd )
    double log10_poisson_cdf_Q_large_lambda ( unsigned int k, double lbd )

from khash cimport *

cdef class P_Score_Upper_Tail:
    """Unit to calculate -log10 poisson_CDF of upper tail and cache
    results in a hashtable.
    """
    cdef:
        kh_int64_t *pscore_table
        
    def __init__ ( self, size_hint = 1 ):
        if size_hint is not None:
            kh_resize_int64( self.pscore_table, size_hint )
            
    def __cinit__( self ):
        self.pscore_table = kh_init_int64()

    def __dealloc__( self ):
        kh_destroy_int64(self.pscore_table)

    cpdef bint check_cache ( self, int x, double l ):
        """Check if the Poisson CDF(x|l) has been calculated and
        cached.
         
        """
        cdef:
            khiter_t k                  # index got from get_init64 function
            long key_value

        key_value = hash( (x, l) )
        k = kh_get_int64( self.pscore_table, key_value )
        return k != self.pscore_table.n_buckets

    cpdef get_pscore ( self, int x, double l ):
        cdef:
            khiter_t k                  # index got in the table; translated from hash key
            int ret = 0
            double val
            long key_value              # hash key
        key_value = hash( (x, l) )
        k = kh_get_int64(self.pscore_table, key_value) # translate hash key to index 
        if k != self.pscore_table.n_buckets: #kh_exist_int64( self.pscore_table, k ):
            return self.pscore_table.vals[k]
        else:
            # calculate and cache
            val = -1 * log10_poisson_cdf ( x, l, 0 )
            k = kh_put_int64( self.pscore_table, key_value, &ret )
            self.pscore_table.keys[ k ] = key_value
            self.pscore_table.vals[ k ] = val
            return val
        
#    def keys( self ):
#        return self.pscore_table.keys
