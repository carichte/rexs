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

class datadict(dict):
    x = None

def norm(vector):
    return np.sqrt((np.array(vector)**2).sum(0))

def lorentzian(x, x0, amp, w):
    return amp/(1 + ((x-x0)/w)**2)

def gaussian(x, x0, amp, w, y0=0):
    return amp*np.exp(-(x-x0)**2/(2*w**2))+y0

def standard_normal(x):
    return np.exp(-x**2/2)/np.sqrt(2*np.pi)

def pvoigt(x, x0, amp, w, y0=0, eta=0.5, feta=None):
    if feta!=None:
        eta = special.erfc(-feta)/2.
    return y0+amp*(eta/(1+((x-x0)/w)**2)+(1-eta)*np.exp(-np.log(2)*((x-x0)/w)**2))

def kasten(x, w):
    return (x>(-w/2.)) * (x<(w/2.))


def lorentzian_filter1d(input, fwhm, axis = -1, output = None,
                      mode = "reflect", cval = 0.0):
    """One-dimensional Lorentzian filter.

    Parameters
    ----------
    %(input)s
    sigma : scalar
        standard deviation for Gaussian kernel
    %(axis)s
    order : {0, 1, 2, 3}, optional
        An order of 0 corresponds to convolution with a Gaussian
        kernel. An order of 1, 2, or 3 corresponds to convolution with
        the first, second or third derivatives of a Gaussian. Higher
        order derivatives are not implemented
    %(output)s
    %(mode)s
    %(cval)s
    """
    fwhm = abs(float(fwhm))
    #sd = fwhm/2.35482
    # make the length of the filter equal to 4 times the standard
    # deviations:
    lw = int(20.0 * fwhm + 0.5)
    weights = [0.0] * (2 * lw + 1)
    weights[lw] = 1.0
    sum = 1.0
    #sd = sd * sd
    # calculate the kernel:
    for ii in range(1, lw + 1):
        #tmp = math.exp(-0.5 * float(ii * ii) / sd)
        #tmp = 1./np.pi * fwhm / (fwhm * fwhm + ii * ii)
        tmp = 1. / (1. + ii * ii / (fwhm * fwhm))
        weights[lw + ii] = tmp
        weights[lw - ii] = tmp
        sum += 2.0 * tmp
    for ii in range(2 * lw + 1):
        weights[ii] /= sum
    return ndimage.correlate1d(input, weights, axis, output, mode, cval, 0)


def plot_fio_files(root, fig=None, log=False, Cols=None):
    """
        Plots all .fio files that contain the string *root*.
        
        Arguments:
        ----------
        
        fig : <class 'matplotlib.figure.Figure'>
            matplotlib figure class to plot data into.
        
        log : bool
            If True, plots in logarithmic scaling in the y axis.
        
        Cols : list
            List of selected Columns to plot.
    """
    if not fig:
        from pylab import figure
        fig=figure()
    FILEDIR=os.path.dirname(root)
    if not FILEDIR: FILEDIR = "."
    MyAxes = {}
    for fioname in os.listdir(FILEDIR):
        if root not in fioname or ".fio" not in fioname: continue
        NextData = FIOdata(FILEDIR + "/" + fioname)
        if not Cols:
            NumCols = len(NextData.data[0])-1
            Cols = range(NumCols)
        else: NumCols = len(Cols)
        for i in range(NumCols):
            i+=1
            if not MyAxes.has_key(i): MyAxes[i] = fig.add_subplot(NumCols,1,i)
            if log: MyAxes[i].semilogy(NextData[:,0], NextData[:,Cols[i-1]])
            else: MyAxes[i].plot(NextData[:,0], NextData[:,Cols[i-1]])
    fig.canvas.draw()
    #fig.canvas.show()
    return MyAxes


def write_header(filename, header):
    """
        Write header (one line) to file
    """
    if not isinstance(header, str):
        raise ValueError("header must be string")
    elif header=="":
        return
    elif not header.startswith("#"):
        header = "#" + header
    if header.endswith("\n") or header.endswith("\r"): pass
    else: header += os.linesep
    myfile = open(filename, "r")
    content = myfile.read()
    myfile.close()
    myfile = open(filename, "w")
    myfile.write(header)
    myfile.write(content)
    myfile.close()


def read_header(filename):
    """
        reads only first line of file.
    """
    myfile = open(filename, "r")
    header = myfile.readline()
    myfile.close()
    return header.strip("\r\n#")


def savedat(fname, data, header="", xcol=None, **kwargs):
    """
        writes data to my standard .dat ascii file type using numpy.savetxt
        and adds a header if given. see numpy.savetxt for further inputs.
        
    """
    if isinstance(data, dict):
        if xcol==None:
            raise ValueError("If dictionary is given, the 1st column (xcol) has to be specified")
        if isinstance(xcol, dict):
            header = " ".join(xcol.keys() + data.keys())
            data = np.vstack(xcol.values() + data.values()).T
        elif isinstance(xcol, str):
            xval = data.pop(xcol)
            header = " ".join([xcol] + data.keys())
            data = np.vstack([xval] + data.values()).T
        
    np.savetxt(fname, data, **kwargs)
    write_header(fname, header)


def loaddat(fname, todict=False):
    """
        open my standard .dat file
    """
    try: 
        data = np.loadtxt(fname, comments=None)
        header = ""
        #print("No header present")
    except:
        data = np.loadtxt(fname, skiprows=1)
        header = read_header(fname)
        cols = header.split()
        if todict and len(cols) == len(data[0]):
            ddata = datadict([(k,data[:,cols.index(k)]) for k in cols])
            ddata.x = cols[0]
            data = ddata
    return data, header


def PolynomialFit(x, y, anchors=None, avgrange=0, order=2, indf=None):
    """
        Returns a fitted polynomial.
        Two scenarios are supported:
            
            1) A set of anchors is given, where the polynomial is being pinned.
               Then the degree equals the amount of anchors given minus 1.
               The y-values for the pins can be averaged over a given range.
               
               if no anchors are given, the following applies:
            
            2) An order of the polynomial is given and the polynomial is
               fitted accordingly. Optionally, a fit range can be given 
               via a masked array ``indf``.
        inputs:
            
            x : array of independent values
            y : array of dependent values
            
            Optional:
                
                anchors  : a set of x-positions where the polynomial
                           shall be pinned to the given curve.
                
                avgrange : float
                           a range in which the given curve will be averaged
                           to calculate the pin position.
                
                
                order : int
                        sets the order of the polynomial for fit.
                        ignored if anchors!=None
                
                indf : boolean array, same length as x and y
                       a mask for x and y to choose a range for fit.
    """
    if anchors==None:
        N = order + 1
        func = lambda x, *v: sum([v[i]*x**(N - 1 - i) for i in reversed(range(N))])
        ind = ~np.isnan(y)
        if indf != None:
            ind *= indf
        start = [1]
        start.extend((N-1)*[0])
        popt, pcov = optimize.curve_fit(func, x[ind], y[ind], p0=start)
        poly = func(x, *popt)
    else:
        if avgrange==0:
            indizes = [abs(xi - x).argmin() for xi in anchors]
        else:
            indizes = [(x > xi - avgrange)*(x < xi + avgrange) for xi in anchors]
        #print sum(indizes).any()
        anchors = np.array([x[ind].mean() for ind in indizes])
        yi      = np.array([y[ind].mean() for ind in indizes])
        
        print anchors, yi
        N = len(indizes)
        v0 =  [0.] * (N - 1) + [1.]
        res = lambda v: sum([v[i]*anchors**(N - 1 - i) for i in range(N)]) - yi
        vopt = optimize.broyden1(res, v0)
        poly = sum([vopt[i]*x**(N - 1 - i) for i in range(N)])
    
    return poly

class TempExpansion(object):
    """
        Class to calculate the anisotropic thermal expansion 
        from given values of the coefficient.
    """
    def __init__(self):
        self.alpha = np.zeros((3,3))
        self.beta  = np.zeros((3,3))
        self.gamma = np.zeros((3,3))
        self.Tref = 298.
    
    def __call__(self, Temperature):
        shrinkage = np.diag((1,1,1)) \
                  + np.array(self.alpha) * (Temperature - self.Tref) \
                  + np.array(self.beta)  * (Temperature - self.Tref)**2 \
                  + np.array(self.gamma) * (Temperature - self.Tref)**3
        return shrinkage

