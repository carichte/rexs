#!/usr/bin/env python
#
#
# this is part of the rexs package written by Carsten Richter
# (carsten.richter@desy.de)
#
import numpy as np
from scipy import optimize
from . import wrap4leastsq

class fitresults(object): pass 

def fitls(x_m, y_m, func, guess, variables="all", power=1, weights=None, 
          fitalg="leastsq", verbose=False, **kwargs):
    """
        
        Fitting a pre-defined function to a set of data using `leastsq' or
        `fmin' from scipy.optimize.
        
        Inputs:
            x_m : numpy.ndarray
                independent values of data points
            
            y_m : numpy.ndarray
                dependent values of data points
            
            func : function
                function to describes the (x_m,y_m) curve. It must accept
                x_m as first argument. The remaining arguments are parameters
                that can be varied to fit `func' onto the data points.
            
            guess : dictionary
                A dictionary containing (parameter, startvalues) items where
                each parameter is one of func's arguments.
            
            variables : list (optional)
                A list of parameter names of `guess' which will be varied to
                fit `func' to the given data (x_m, y_m).
                default: guess.keys() (all parameters)
            
            power : float (optional)
                The differences between data points and fit function can be 
                raised to this power to modify the influence which small or 
                large devieations have.
                default: 1
            
            weights : str or numpy.ndarray (optional)
                Defines how the different data points shall be weighted 
                during fit. I can either be an array containing the weight for
                each of the data points or 'statistical' to weight by the
                inverse square root of the data values `y_m' like in poisson
                statistics or, last, 'relative' to weight by the inverse of
                the data values if the relative deviation of func and data 
                is to be minimized.
                Default: no weighting
            
            fitalg : str (optional)
                Can be either 'leastsq' or 'simplex' to make use of the 
                scipy.optimize functions `leastsq' or `fmin', respectively.
                Default: 'leastsq'
        
        Returns:
            A `fitresults' object containing the following attributes:
            
            - fitresult.popt : dictionary
                A dictionary like `guess' but containing the values of the 
                parameters after optimization
            - fitresult.err : float
                The value of the used residuals function in the fount optimum
            - fitresult.stddev : dictionary
                A dictionary like `guess' containing the standard deviation 
                of each parameter in case `fitalg="leastsq"' was used.
            - fitresult.Rval : float
                The R-value describing the goodness of fit as is used in 
                X-ray diffraction.
            - fitresult.ysim : numpy.ndarray
                The evaluated values of func at `x_m' in the optimum
            - fitresult.yguess : numpy.ndarray
                The evaluated values of func at `x_m' for the guess
            - fitresult.func : function
                The fitfunction `func' for the found solution using the 
                optimized parameters and only accepting `x_m' as an input.
            
        
        Example: See ../examples/fitls/fitexample.py
        
    """
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
            return res
    elif weights=="relative":
        ind = y_m>0
        def residuals(v):
            return (func(x_m[ind], **v) / y_m[ind] - 1)**power
    elif weights=="Rval":
        def residuals(v):
            return (func(x_m,**v) - y_m)/abs(y_m).sum()
    elif hasattr(weights, "__iter__")\
     and hasattr(weights, "__len__")\
     and len(weights)==len(x_m):
        def residuals(v):
            return (func(x_m, **v) - y_m)**power * weights
    elif weights is None:
        def residuals(v):
            return (func(x_m, **v) - y_m)**power
    else: 
        def residuals(v):
            return (func(x_m, **v) - y_m)**power
        print("Invalid input for ``weights``. No weighting implemented.")
    
    
    def residualswrap(v):
        if verbose:
            err = residuals(v)
            print("%g"%(err**2).sum())/len(x_m)
            return err
        else:
            return residuals(v)
    sumofsquares = lambda v: (residualswrap(v)**2).sum()
    
    
    fitted_param = guess.copy()
    if variables=="all":
        variables = list(guess)
    if fitalg=="leastsq":
        fitfunction, startvalues = wrap4leastsq.wrap_for_fit(residualswrap, guess, 
                                                      variables, unpack=False)
        output = optimize.leastsq(fitfunction, startvalues, full_output=True,**kwargs)
        #print output
        if output[1] is None: 
            stddev = [np.inf for i in range(len(variables))]
        else: 
            stddev = [np.sqrt(var) for var in output[1].diagonal()] # Kovarianzmatrix
    elif fitalg=="simplex":
        fitfunction, startvalues = wrap4leastsq.wrap_for_fit(sumofsquares, 
                                              guess, variables, unpack=False)
        if not "maxfun" in kwargs:
            kwargs["maxfun"] =  500*len(startvalues)
        if not "maxiter" in kwargs:
            kwargs["maxiter"] = 500*len(startvalues)
        output = optimize.fmin(fitfunction, startvalues, full_output=True, **kwargs) 
        
        stddev = [np.inf for i in range(len(variables))]
    if len(variables)>1: 
        fitted_param.update(dict(zip(variables, output[0]))) # minimize
        stddev_dict = dict(zip(variables, stddev))
    else:
        fitted_param.update(dict([(variables[0], output[0])])) # minimize
        stddev_dict = dict([(variables[0], stddev)])
    
    result = fitresults()
    result.yguess = func(x_m, **guess)
    y_s = func(x_m, **fitted_param)
    result.popt = fitted_param
    result.err = sumofsquares(fitted_param)/len(x_m)
    result.stddev = stddev_dict
    result.Rval = abs(y_m - y_s).sum()/abs(y_m).sum()
    result.ysim = y_s
    result.func = lambda x: func(x, **fitted_param)
    if verbose:
        print("Results of %s fit:\n"    
              " sum of residuals: %g\n"
              " Rval: %g"%(fitalg, result.err, result.Rval))
    return result
    


