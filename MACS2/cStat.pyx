# Time-stamp: <2012-02-07 20:08:05 Tao Liu>

"""Module Description

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

from collections import Counter

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
def histogram ( vl, breaks=None, minv=None, maxv=None, binsize=None):
    """Return histogram statistics.

    Parameters:

    vl: 2D numpy.array as [ [value, length], [value, length], ...]
    
    breaks: if breaks is not None and a valid integar, split [min,max]
    of values in vl into number of equal sized bins. Otherwise, no
    binning is involved.

    Return Value:
    Counter object


    when breaks is not None, key values in Counter is the start points
    of each bin.
    
    """
    assert breaks == None or isinstance(breaks,int)
    
    ret = Counter()

    if breaks == None and binsize == None:
        for (v,l) in vl:
            ret[v] += int(l)
    else:
        if maxv == None:
            maxv = vl[:,0].max()
        if minv == None:
            minv = vl[:,0].min()
        if binsize == None:
            binsize = (maxv-minv)/breaks
        for (v,l) in vl:
            k = (v - minv)//binsize*binsize + minv
            #print k
            ret[ k ] += int(l)

    return ret

def histogram2D ( md ):
    """Return histogram statistics.

    Parameters:

    vl: 2D numpy.array as [ [value, length], [value, length], ...]
    
    breaks: if breaks is not None and a valid integar, split [min,max]
    of values in vl into number of equal sized bins. Otherwise, no
    binning is involved.

    Return Value:
    Counter object


    when breaks is not None, key values in Counter is the start points
    of each bin.
    
    """
    ret = Counter()

    for (m, d, l) in md:
        ret[ (m,d) ] += int(l)

    return ret

# ------------------------------------
# Classes
# ------------------------------------
