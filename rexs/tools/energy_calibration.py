#!/usr/bin/env python
#
#
# this is part of the rexs package written by Carsten Richter
# (carsten.richter@desy.de)
#
import os

import numpy as np
from scipy import optimize, interpolate

def energy_calibration(master, slave, col_ref, drift=False,
                              splorder=3, dE=0, indf=None):
    """
        Function to perform x-axis correction for measurements containing
        data measured from a non changing reference. Frequently used in
        absorption spectroscopy where x equals energy.
        The correction involves a shift of the data contained in slave
        to match the reference curves contained in in column `col_ref`.
        The data contained in master remains untouched.
        
        Inputs:
            - master : numpy.array
                2D data array containing x-values in 1st column and y values
                including the one with reference data in the remaining
                columns.
            
            - slave : numpy.array
                similar as master but y-data will be shifted to match the
                reference columnes
            
            - col_ref : integer or 2-tuple
                defines which columns are containing the reference data in
                the (master, slave) data array, respectively.
            
            - drift : bool
                if True, a linear x-dependence will be considered in
                addition to the constant shift in x.
            
            - splorder : int in [1,2,3]
                defines the polynomial order which is being used to
                interpolate the slave data curves for
            
            - dE : float
                the initial guess for the difference in x-values between
                slave and master.
            
            - indf - numpy.array as bool
                1D mask array to select values that shall be
                considered in the offset correction.
                must have same amount of lines as does master and err
    """
    if hasattr(col_ref, "__iter__"):
        refm = col_ref[0]
        refs = col_ref[1]
    else:
        refm = col_ref
        refs = col_ref
    
    if indf is None:
        indf = 1
    else:
        assert master.shape[0]==indf.shape[0],\
            "lengths of master slave and indf have to agree"
    
    
    slave = slave.copy()
    ind = (master[:,0] > slave[:,0].min()) * (master[:,0] < slave[:,0].max())
    ind *= indf
    if not ind.any():
        return slave,[]
    
    slavefunc = interpolate.UnivariateSpline(slave[:,0], slave[:,refs],
                                                       k=splorder, s=0)
    #slavefunc = interpolate.interp1d(slave[:,0], slave[:,refs],
    #                 kind=splorder, bounds_error=False, fill_value=np.nan)
    def newslave(offset, scale, drift):
        return scale*slavefunc(master[ind,0]*drift + offset)
    def residuals(*v):
        if hasattr(v[0], "__iter__"): 
            v = np.append(v[0], v[1:]).ravel()
        #return (np.diff(newslave(*v)/master[ind,refm]-1))
        #return ((newslave(*v)/master[ind,col_ref]-1))
        return (newslave(*v)-master[ind,refm])
    #fp = optimize.fmin(residuals, [0], args=(1,))
    fp = optimize.leastsq(residuals, [dE,1.], args=(1.,))[0]
    #if hasattr(fp[0], "__iter__"): fp = np.append(fp[0], fp[1:]).ravel()
    #print fp
    newx = (slave[:,0] - fp[0])
    if drift:
        fp = np.append(fp,1.)
        #fp = optimize.fmin(residuals, [0,1])
        fp = optimize.leastsq(residuals, fp)[0]
        #if hasattr(fp[0], "__iter__"): fp = np.append(fp[0], fp[1:]).ravel()
        #print fp
        newx = (slave[:,0] - fp[0])/fp[2]
    slavefunc = interpolate.interp1d(newx, slave[:,1:].T, kind=splorder,
                                     bounds_error=False, fill_value=np.nan)
    newslave = np.vstack((slave[:,0], slavefunc(slave[:,0]))).T
    #print fp
    fp = dict(zip(["offset", "scale", "drift"],fp))
    out = dict({"popt":fp, "func":slavefunc})
    return newslave, out


