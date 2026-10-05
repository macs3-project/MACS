# cython: language_level=3

# Py_REFCNT from Python.h. Cython 3.1 and later declare it in
# cpython/ref.pxd; Cython 3.0 does not, so it is declared here.
cdef extern from "Python.h":
    Py_ssize_t Py_REFCNT(object o)
