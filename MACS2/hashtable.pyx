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

    cpdef get_item(self, int64_t key):
        """Given a key, return a value.
        
        """
        cdef khiter_t k
        k = kh_get_int64(self.table, key)
        if k != self.table.n_buckets:
            return self.table.vals[k]
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
        self.table = kh_init_float64()

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

# cdef class PyObjectHashTable:

#     cdef:
#         kh_pymap_t *table

#     def __init__(self, size_hint=1):
#         self.table = kh_init_pymap()
#         kh_resize_pymap(self.table, size_hint)

#     def __dealloc__(self):
#         if self.table is not NULL:
#             self.destroy()

#     cpdef destroy(self):
#         kh_destroy_pymap(self.table)
#         self.table = NULL

#     cpdef get_item(self, object val):
#         cdef khiter_t k
#         k = kh_get_pymap(self.table, <PyObject*>val)
#         if k != self.table.n_buckets:
#             return self.table.vals[k]
#         else:
#             raise KeyError(val)

#     def get_iter_test(self, object key, Py_ssize_t iterations):
#         cdef Py_ssize_t i, val
#         for i in range(iterations):
#             k = kh_get_pymap(self.table, <PyObject*>key)
#             if k != self.table.n_buckets:
#                 val = self.table.vals[k]

#     cpdef set_item(self, object key, Py_ssize_t val):
#         cdef:
#             khiter_t k
#             int ret = 0
#             char* buf

#         k = kh_put_pymap(self.table, <PyObject*>key, &ret)
#         # self.table.keys[k] = key
#         if kh_exist_pymap(self.table, k):
#             self.table.vals[k] = val
#         else:
#             raise KeyError(key)

#     def map_locations(self, ndarray[object] values):
#         cdef:
#             Py_ssize_t i, n = len(values)
#             int ret = 0
#             object val
#             khiter_t k

#         for i in range(n):
#             val = values[i]
#             k = kh_put_pymap(self.table, <PyObject*>val, &ret)
#             self.table.vals[k] = i

#     def lookup(self, ndarray[object] values):
#         cdef:
#             Py_ssize_t i, n = len(values)
#             int ret = 0
#             object val
#             khiter_t k
#             ndarray[int32_t] locs = np.empty(n, dtype='i4')

#         for i in range(n):
#             val = values[i]
#             k = kh_get_pymap(self.table, <PyObject*>val)
#             if k != self.table.n_buckets:
#                 locs[i] = self.table.vals[k]
#             else:
#                 locs[i] = -1

#         return locs

#     def lookup2(self, ndarray[object] values):
#         cdef:
#             Py_ssize_t i, n = len(values)
#             int ret = 0
#             object val
#             khiter_t k
#             long hval
#             ndarray[int32_t] locs = np.empty(n, dtype='i4')

#         # for i in range(n):
#         #     val = values[i]
#             # hval = PyObject_Hash(val)
#             # k = kh_get_pymap(self.table, <PyObject*>val)

#         return locs

#     def unique(self, ndarray[object] values):
#         cdef:
#             Py_ssize_t i, n = len(values)
#             Py_ssize_t idx, count = 0
#             int ret = 0
#             object val
#             khiter_t k
#             list uniques = []
#             bint seen_na = 0

#         for i in range(n):
#             val = values[i]

#             if not _checknan(val):
#                 k = kh_get_pymap(self.table, <PyObject*>val)
#                 if k == self.table.n_buckets:
#                     k = kh_put_pymap(self.table, <PyObject*>val, &ret)
#                     uniques.append(val)
#             elif not seen_na:
#                 seen_na = 1
#                 uniques.append(ONAN)

#         return uniques

#     cpdef get_labels(self, ndarray[object] values, list uniques,
#                      Py_ssize_t count_prior, int32_t na_sentinel):
#         cdef:
#             Py_ssize_t i, n = len(values)
#             ndarray[int32_t] labels
#             ndarray[int32_t] counts
#             Py_ssize_t idx, count = count_prior
#             int ret = 0
#             object val
#             khiter_t k

#         labels = np.empty(n, dtype=np.int32)
#         counts = np.empty(count_prior + n, dtype=np.int32)

#         for i in range(n):
#             val = values[i]

#             if val != val or val is None:
#                 labels[i] = na_sentinel
#                 continue

#             k = kh_get_pymap(self.table, <PyObject*>val)
#             if k != self.table.n_buckets:
#                 idx = self.table.vals[k]
#                 labels[i] = idx
#                 counts[idx] = counts[idx] + 1
#             else:
#                 k = kh_put_pymap(self.table, <PyObject*>val, &ret)
#                 self.table.vals[k] = count
#                 uniques.append(val)
#                 labels[i] = count
#                 counts[count] = 1
#                 count += 1

#         return labels, counts[:count].copy()

#     # def unique(self, ndarray[object] values, list uniques):
#     #     cdef:
#     #         Py_ssize_t i, n = len(values)
#     #         Py_ssize_t idx, count = 0
#     #         int ret
#     #         object val
#     #         khiter_t k

#     #     for i in range(n):
#     #         val = values[i]
#     #         k = kh_get_pymap(self.table, <PyObject*>val)
#     #         if k == self.table.n_buckets:
#     #             k = kh_put_pymap(self.table, <PyObject*>val, &ret)
#     #             uniques.append(val)
#     #             count += 1

