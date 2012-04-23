# Time-stamp: <2012-04-23 04:53:32 Tao Liu>

"""Module Description: Fast/smaller array using C code

Copyright (c) 2012 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

cdef extern from *:
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"
    
cdef extern from "stdlib.h":
    void free(void* ptr)
    void* malloc(size_t size)
    void* realloc(void *ARRAY, size_t SIZE)
    void qsort (void *ARRAY, size_t COUNT, size_t SIZE, int (*COMPARE)(const_void *, const_void *)) nogil


# ------------------------------------
# python modules
# ------------------------------------

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

cdef int compare_int(const_void *a, const_void *b):
    return (<int *>a)[0] - (<int *>b)[0]

cdef int compare_float(const_void *a, const_void *b):
    return 1 if (<float *>a)[0] > (<float *>b)[0] else 0

# ------------------------------------
# Classes
# ------------------------------------

cdef class IntArray:
    cdef int * d
    cdef long ia
    cdef long ii
    cdef long m
    
    def __init__ (self, long l):
        self.ia = 0                      # for index to add data/last index of valid data
        self.ii = 0                      # for index to get data
        self.m = l                       # maximum index
        cdef int * data = <int *>malloc( l*sizeof(int) )
        self.d = data
        
    def __dealloc__ (self):
        free( <void *> self.d )
        self.d = NULL
        
    cpdef put (self, int a):
        assert self.ia < self.m, "Out of boundary"
        #if self.ia == self.m:
        #    raise IndexError
        self.d[self.ia] = a
        self.ia += 1

    cpdef put_i (self, long i, int a):
        assert i < self.m, "Out of boundary!"
        self.d[i] = a

    cpdef get (self, long i):
        assert i < self.m, "Out of boundary!"
        return self.d[i]

    cpdef sort (self):
        qsort(self.d, self.ia, sizeof(int) , compare_int)

    cpdef resize (self, long new_size):
        self.m = new_size
        self.d = <int *>realloc( self.d, new_size*sizeof(int) )

    def __iter__ (self):
        self.ii = -1
        return self

    def __next__ (self):
        if self.ii == self.m:
            raise StopIteration
        self.ii = self.ii + 1
        return self.d[self.ii]

cdef class FloatArray:
    cdef float * d
    cdef long ia
    cdef long ii
    cdef long m
    
    def __init__ (self, long l):
        self.ia = 0                      # for index to add data
        self.ii = 0                      # for index to get data
        self.m = l                      # maximum index
        cdef float * data = <float *>malloc( l*sizeof(float) )
        self.d = data
        
    def __dealloc__ (self):
        free( <void *> self.d )
        self.d = NULL
        
    cpdef put (self, float a):
        assert self.ia < self.m, "Out of boundary"
        #if self.ia == self.m:
        #    raise IndexError
        self.d[self.ia] = a
        self.ia += 1

    cpdef put_i (self, long i, float a):
        assert i < self.m, "Out of boundary!"
        self.d[i] = a

    cpdef get (self, long i):
        assert i < self.m, "Out of boundary!"
        return self.d[i]

    cpdef sort (self):
        qsort(self.d, self.m, sizeof(float) , compare_float)

    cpdef resize (self, long new_size):
        self.m = new_size
        self.d = <float *>realloc( self.d, new_size*sizeof(float) )

    def __iter__ (self):
        self.ii = -1
        return self

    def __next__ (self):
        if self.ii == self.m:
            raise StopIteration
        self.ii = self.ii + 1
        return self.d[self.ii]

