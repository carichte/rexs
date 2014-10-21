"""
Filename: rebin.py
    (1) rebin1d(x,y,newx)
"""
import os
import numpy as nm
import ctypes as ct

_librebin = nm.ctypeslib.load_library('librebin', os.path.dirname(__file__))


_librebin.rebin.argtypes = [nm.ctypeslib.ndpointer(dtype = nm.double),
                            nm.ctypeslib.ndpointer(dtype = nm.double),
                            nm.ctypeslib.ndpointer(dtype = nm.double),
                            nm.ctypeslib.ndpointer(dtype = nm.double),
                            ct.c_int,
                            ct.c_int
                            ]

_librebin.rebin_from_centers.argtypes = [
                            nm.ctypeslib.ndpointer(dtype = nm.double),
                            nm.ctypeslib.ndpointer(dtype = nm.double),
                            nm.ctypeslib.ndpointer(dtype = nm.double),
                            nm.ctypeslib.ndpointer(dtype = nm.double),
                            ct.c_int,
                            ct.c_int
                            ]

def rebin1d(x, y, newx, average=False):
    """
        Computes the weighted rebinning of density dirstribution data taking 
        into account partial overlap of the original and new intervals.
        
        Usually the borders of the bins `x' are given and the new values are
        calculated for the given borders `newx' of the new bins. In this case,
        the input array x is by one element longer than y.
        
        However, this function also supports the case when the input `x' 
        contains the centers of the bins and, thus, is as long as `y'.
        Then, the new values are also calculated for the given centers `newx'
        of the new bins.
        
        ARGUMENT(S):
            x - 1darray, float
                the centers or borders of the original bins
            y - 1darray, float
                the values of the original bins
            newx - 1darray, float
                   the centers or borders of the new bins
    
        RESULT(S):
            newx - 1darray, float
                   the values of the new bins
    
    """
    x = nm.array(x, dtype=nm.double, ndmin=1)
    y = nm.array(y, dtype=nm.double, ndmin=1)
    newx = nm.array(newx, dtype=nm.double, ndmin=1)
    assert x.ndim == 1 and y.ndim==1 and newx.ndim == 1
    assert (nm.diff(x)>0).all() and (nm.diff(newx)>0).all(),\
        "Both x and newx must be ordered."
    if len(x) == len(y)+1:
        newy = nm.empty(len(newx)-1, dtype=nm.double)
        _librebin.rebin(x,y,newx,newy, len(y), len(newy))
        if average:
            num = nm.empty(len(newx)-1, dtype=nm.double)
            _librebin.rebin(x, nm.ones_like(y), newx, num, len(y), len(num))
            newy/=num
    elif len(x) == len(y):
        newy = nm.empty_like(newx, dtype=nm.double)
        _librebin.rebin_from_centers(x,y,newx,newy, len(y), len(newy))
        if average:
            num = nm.empty_like(newx, dtype=nm.double)
            _librebin.rebin_from_centers(x, nm.ones_like(y), newx, num, 
                                            len(y), len(num))
            newy/=num
    else:
        raise ValueError("Invalid lengths of x and y")
    return newy
