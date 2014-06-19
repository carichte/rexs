#
#
# this is part of the 'evaluationtools' set written by Carsten Richter
# (carsten.richter@desy.de)
#

import numpy as np
import os
from scipy import ndimage



def BackgroundXRD(pattern, span=20, ftol=1e-4, maxiter=2500, nsigma=1, 
                  verbose=False):
    """
        Function to remove Background from a powder XRD pattern
        
        Input:
            pattern : 1D numpy.ndarray
                intensity values
            span : int
                half width of smoothing window in channels
            ftol : float
                goal tolerance in relative deviation
            maxiter : int
                maximum number of iterations
            nsigma : float
                tolerance in terms of standard deviations for distinction of 
                background and signal (modification)
            verbose : bool
                print extra output?
        
        The algorithm is a small modification of the one described by 
        Sergio Bruckner in:
            J. Appl. Cryst. (2000). 33, 977-979 
            http://dx.doi.org/10.1107/S0021889800003617
    """
    pattern = pattern.copy()
    
    pmean = pattern.mean()
    pconst = pmean + 2*(pmean - pattern.min())
    ind = pattern > pconst
    pattern[ind] = pconst
    
    pleft = pattern[:span].min()
    pright = pattern[-span:].min()
    pattern = np.append(np.ones(span)*pleft, 
                        np.append(pattern, np.ones(span)*pright))
    err = np.inf
    derr= 1.
    i=0
    while derr>ftol and i<maxiter:
        #if not ((i+1)%100):
        #    print "median"
        #    pattern =  ndimage.median_filter(pattern, 3)
        smoothpattern = ndimage.uniform_filter1d(pattern, size=(2*span+1),
                                                 mode="reflect")
        newerr=((smoothpattern - pattern)**2).sum()
        ind = (pattern + nsigma*np.sqrt(pattern))>smoothpattern
        pattern[ind]=smoothpattern[ind]
        
        derr = abs(newerr/err-1.)
        err = newerr
        i+=1
        #print i, derr
    if verbose:
        print "Number of Iterations: %i"%i
        print "Final Relative Error: %g"%derr
    return pattern[span:-span]
