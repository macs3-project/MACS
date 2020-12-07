# cython: language_level=3
# cython: profile=True
# cython: linetrace=True
# Time-stamp: <2020-12-06 21:32:47 Tao Liu>
import numpy as np
cimport numpy as np
from numpy cimport ndarray
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t
import multiprocessing as mp
from MACS3.Signal.Prob import poisson_cdf

def cal_pscore_mp ( int32_t N, ndarray pos, ndarray treat, ndarray ctrl, dict p_dict ):
    #cdef:
    #    int32_t n_queries, n_mp_range, i
    #    int32_t s, e, p0
    #    list args
    #    
    n_queries = pos.shape[ 0 ]
    # decide how to split the tasks to N_mp slices
    n_mp_range = n_queries // N + 1
    # [0*n_mp_range .. 1*n_mp_range] [1*n_mp_range .. 2*n_mp_range] ...
    # start to pool
    P = mp.Pool( N )
    args = []
    # prepare args
    for i in range( N ):
        s = i * n_mp_range
        e = (i + 1) * n_mp_range
        if e > n_queries:
            e = n_queries
        if s == 0:
            p0 = 0
        else:
            p0 = pos[ s-1 ]
        args.append( ( p_dict, p0, pos, treat, ctrl, s, e ) )
    mapresults = P.map_async( __update_pscore, args )
    P.close()
    P.join()
    results = mapresults.get() # for each result in results,
                                   # result[0] is partial pscore_dict,
                                   # result[1] is partial pscore_stat
                                   # now merge partial results to full
                                   # results: 1. global pscore_dict,
                                   # 2. own pvalue_stat
    return results
                                       
cpdef __update_pscore ( tuple x ):
    cdef:
        dict p_dict
        int32_t p0
        ndarray p
        ndarray t
        ndarray c
        int32_t s
        int32_t e
        float32_t v
        int32_t *p_ptr
        float32_t *t_ptr
        float32_t *c_ptr
        dict pscore_stat
        dict extra_p_dict
        #int i
        #tuple tu

    extra_p_dict = {}
    pscore_stat = {}
    ( p_dict, p0, p, t, c, s, e ) = x
    pscore_stat = {}
    p_ptr = <int32_t *> p.data
    t_ptr = <float32_t *> t.data
    c_ptr = <float32_t *> c.data
    p_ptr += s
    t_ptr += s
    c_ptr += s
    
    pre_p = p0
    for i in range( e - s ):
        tmp_tuple = ( <int32_t> ( t_ptr[0] ), c_ptr[0] )
        if tmp_tuple in p_dict:
            v = p_dict[ tmp_tuple ]
        else:
            # calculate and cache
            v = -1.0 * poisson_cdf ( tmp_tuple[0], tmp_tuple[1], False, True )
            p_dict[ tmp_tuple ] = v
            extra_p_dict[ tmp_tuple ] = v 
        l = p_ptr[0] - pre_p
        if v in pscore_stat:
            pscore_stat[ v ] += l
        else:
            pscore_stat[ v ] = l
        pre_p = p_ptr[0]
        p_ptr += 1
        t_ptr += 1
        c_ptr += 1
    return ( extra_p_dict, pscore_stat )
