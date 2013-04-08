"""
Modified from Pandas: https://github.com/pydata/pandas/blob/master/pandas/src/hashtable.pyx
"""

from khash cimport *
from numpy cimport *
import numpy as np
cimport numpy as np

cdef inline bint _checknan(object val):
    return not np.PyArray_Check(val) and val != val

# def list_to_object_array(list obj):
#     """
#     Convert list to object ndarray. Seriously can't believe I had to write this
#     function
#     """
#     cdef:
#         Py_ssize_t i, n
#         ndarray[object] arr

#     n = len(obj)
#     arr = np.empty(n, dtype=object)

#     for i from 0 <= i < n:
#         arr[i] = obj[i]

#     return arr


cdef class Int64HashTable:
    """A hashtable taking 64bit integer as key and 64bit float as
    value.

    """
    cdef:
        kh_int64_t *table

    def __init__(self, size_hint=1):
        if size_hint is not None:
            kh_resize_int64(self.table, size_hint)

    def __cinit__(self):
        self.table = kh_init_int64()

    def __dealloc__(self):
        kh_destroy_int64(self.table)

    cpdef bint has_key(self, int64_t key):
        """Check if a given key is valid in hashtable.
        
        """
        cdef khiter_t k
        k = kh_get_int64(self.table, key)
        return k != self.table.n_buckets

    def __getitem__(self, int64_t key):
        """Given a key, return a value.
        
        """
        cdef khiter_t k
        k = kh_get_int64(self.table, key)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(key)

    cpdef get_item(self, int64_t key):
        """Given a key, return a value.
        
        """
        cdef khiter_t k
        k = kh_get_int64(self.table, key)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(key)

    def __setitem__ ( self, int64_t key, float64_t val):
        """Put a key-value pair to hashtable.
        
        """
        cdef:
            khiter_t k
            int ret = 0

        k = kh_put_int64(self.table, key, &ret)
        self.table.keys[k] = key
        if kh_exist_int64(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)    

    cpdef set_item(self, int64_t key, float64_t val):
        """Put a key-value pair to hashtable.
        
        """
        cdef:
            khiter_t k
            int ret = 0

        k = kh_put_int64(self.table, key, &ret)
        self.table.keys[k] = key
        if kh_exist_int64(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)

    def map(self, ndarray[int64_t] keys, ndarray[float64_t] values):
        """Take an array of keys in int64, and an array of values in
        float64 (of the same length), create key-value pairs in
        hashtable.
        
        """
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t key
            khiter_t k

        for i in range(n):
            key = keys[i]
            k = kh_put_int64(self.table, key, &ret)
            self.table.vals[k] = <float64_t> values[i]

    def map_locations(self, ndarray[int64_t] values):
        """Take a list of keys, set the value of each key according to
        the index in input array.
        
        """
        cdef:
            Py_ssize_t n = len(values)
            float64_t i
            int ret = 0
            int64_t key
            khiter_t k

        for i in range(n):
            key = values[i]
            k = kh_put_int64(self.table, key, &ret)
            self.table.vals[k] = i

    def lookup(self, ndarray[int64_t] values):
        """Take a list of keys, return their values.
        
        """
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t key
            khiter_t k
            ndarray[float64_t] locs = np.empty(n, dtype='float64')

        for i in range(n):
            key = values[i]
            k = kh_get_int64(self.table, key)
            if k != self.table.n_buckets:
                locs[i] = self.table.vals[k]
            else:
                locs[i] = -1

        return locs

cdef class Float64HashTable:
    """A hashtable taking 64bit float as key and 64bit float as
    value.

    """
    cdef:
        kh_float64_t *table

    def __init__(self, size_hint=1):
        if size_hint is not None:
            kh_resize_float64(self.table, size_hint)

    def __cinit__(self):
        self.table = kh_init_float64() # initiate key table with float64 as keys

    def __dealloc__(self):
        kh_destroy_float64(self.table)

    cpdef bint has_key(self, float64_t key):
        """Check if a given key is valid in hashtable.
        
        """
        cdef khiter_t k
        k = kh_get_float64(self.table, key)
        return k != self.table.n_buckets

    cpdef get_item(self, float64_t key):
        """Given a key, return a value.
        
        """
        cdef khiter_t k
        k = kh_get_float64(self.table, key)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(key)

    def __getitem__ ( self, float64_t key ):
        cdef khiter_t k
        k = kh_get_float64(self.table, key)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(key)        

    cpdef set_item(self, float64_t key, float64_t val):
        """Put a key-value pair to hashtable.
        
        """
        cdef:
            khiter_t k
            int ret = 0

        k = kh_put_float64(self.table, key, &ret)
        self.table.keys[k] = key
        if kh_exist_float64(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)

    def __setitem__ (self, float64_t key, float64_t val):
        """Put a key-value pair to hashtable.
        
        """
        cdef:
            khiter_t k
            int ret = 0

        k = kh_put_float64(self.table, key, &ret)
        self.table.keys[k] = key
        if kh_exist_float64(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)

    def map(self, ndarray[float64_t] keys, ndarray[float64_t] values):
        """Take an array of keys in float64, and an array of values in
        float64 (of the same length), create key-value pairs in
        hashtable.
        
        """
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            float64_t key
            khiter_t k

        for i in range(n):
            key = keys[i]
            k = kh_put_float64(self.table, key, &ret)
            self.table.vals[k] = <float64_t> values[i]

    def map_locations(self, ndarray[float64_t] values):
        """Take a list of keys, set the value of each key according to
        the index in input array.
        
        """
        cdef:
            Py_ssize_t n = len(values)
            float64_t i
            int ret = 0
            float64_t key
            khiter_t k

        for i in range(n):
            key = values[i]
            k = kh_put_float64(self.table, key, &ret)
            self.table.vals[k] = i

    def lookup(self, ndarray[float64_t] values):
        """Take a list of keys, return their values.
        
        """
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            float64_t key
            khiter_t k
            ndarray[float64_t] locs = np.empty(n, dtype='float64')

        for i in range(n):
            key = values[i]
            k = kh_get_float64(self.table, key)
            if k != self.table.n_buckets:
                locs[i] = self.table.vals[k]
            else:
                locs[i] = -1

        return locs



ONAN = np.nan

