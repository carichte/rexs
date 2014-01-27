#!/usr/bin/env python
#
#
# this is part of the 'evaluationtools' set written by Carsten Richter
# (carsten.richter@desy.de)
#
import numpy as np
import os
from scipy import ndimage, optimize, special
import time
import wrap4leastsq

class expando(object): pass 

def fitls(x_m, y_m, func, guess, variables, power=1, weights=None, fitalg="leastsq"):
    """
        Fits function `func` to measured data (x_m, y_m) using scipy`s `leastsq` or `fmin`
        
        'guess' is a dictionary of all parameters 'func' takes including starting values.
        'variables' is a list of strings containing the parameters which shall be varied.
        
        Returns a (resulting function, resulting parameters, sum of squares) tuple.
        
        weights can be 'statistical', an array of length len(x_m) or None
    """
    #print hasattr(weights, "__iter__") , hasattr(weights, "__len__") , len(weights)==len(x_m)
    if not isinstance(x_m, np.ndarray):
        x_m = np.array(x_m)
    if not isinstance(y_m, np.ndarray):
        y_m = np.array(y_m)
    if weights=="statistical":
        ind = y_m>0
        counts = y_m/y_m[ind].min()
        weights = np.ones_like(counts)
        weights[ind] = 1./np.sqrt(counts[ind])
        def residuals(v):
            res = (func(x_m, **v) - y_m)**power * weights
            print (res**2).sum()
            return res
    elif weights=="relative":
        ind = y_m>0
        def residuals(v):
            return (func(x_m[ind], **v) / y_m[ind] - 1)**power
    elif hasattr(weights, "__iter__")\
     and hasattr(weights, "__len__")\
     and len(weights)==len(x_m):
        def residuals(v):
            return (func(x_m, **v) - y_m)**power * weights
    elif weights==None:
        def residuals(v):
            return (func(x_m, **v) - y_m)**power
    else: 
        def residuals(v):
            return (func(x_m, **v) - y_m)**power
        print("Bad input for ``weights``. No weighting implemented.")
    
    sumofsquares = lambda v: (residuals(v)**2).sum()
    
    fitted_param = guess.copy()
    if fitalg=="leastsq":
        fitfunction, startvalues = wrap4leastsq.wrap_for_fit(residuals, guess, variables, unpack=False)
        output = optimize.leastsq(fitfunction, startvalues, full_output=True)
        #print output
        if output[1]==None: stddev = [np.inf for i in range(len(variables))]
        else: stddev = [np.sqrt(var) for var in output[1].diagonal()] # Kovarianzmatrix
    elif fitalg=="simplex":
        fitfunction, startvalues = wrap4leastsq.wrap_for_fit(sumofsquares, guess, variables, unpack=False)
        output = optimize.fmin(fitfunction, startvalues, full_output=True, maxfun=1000*len(startvalues), maxiter=1000*len(startvalues))
        stddev = [np.inf for i in range(len(variables))]
    if len(variables)>1: 
        fitted_param.update(dict(zip(variables, output[0]))) # minimize
        stddev_dict = dict(zip(variables, stddev))
    else:
        fitted_param.update(dict([(variables[0], output[0])])) # minimize
        stddev_dict = dict([(variables[0], stddev)])
    y_s = func(x_m, **fitted_param)
    outdict = {}
    result = expando()
    result.popt = fitted_param
    result.err = (residuals(fitted_param)**2).sum()/len(x_m)
    result.stddev = stddev_dict
    result.Rval = abs(y_m - y_s).sum()/abs(y_m).sum()
    result.ysim = y_s
    result.yguess = func(x_m, **guess)
    result.func = lambda x: func(x, **fitted_param)
    return result
    

