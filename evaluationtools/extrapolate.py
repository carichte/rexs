import numpy as np
from scipy import ndimage
from fitls import fitls


def exponential(x, y, ind, tau=1):
    """
        simple axtrapolation assuming exponential decay
        
        ind is a mask for the data to use.
    """
    x = x[ind]
    y = y[ind]
    #ysmooth = ndimage.uniform_filter1d(y, int(len(y)/10.))
    # get some starting parameters
    print int(1./3. * len(y)), y[10]
    d1 = np.diff(y[:(int(1./3. * len(y)))]).mean()
    d2 = np.diff(y[(int(2./3 * len(y))):]).mean()
    if d2>d1:
        if abs(d2)>abs(d1):
            a0 = y.max() - y.min()
            t0 = x[-1] - x[0]
        else:
            a0 = y.max() - y.min()
            t0 = -(x[-1] - x[0])
    else:
        if abs(d2)>abs(d1):
            a0 = y.min() - y.max()
            t0 = x[-1] - x[0]
        else:
            a0 = y.min() - y.max()
            t0 = -(x[-1] - x[0])
            
    fp0 = {"c":x[0], "t":t0, "a":a0, "y0":np.median(y)}
    func = lambda x, a, t, c, y0 : a * np.exp((x - c)/t) + y0
    fp, er, stddev = fitls(x, y, func, fp0, ["a", "y0"],
                         fitalg="leastsq", weights=None)
    #fp, er, stddev = fitls(x, y, func, fp, ["a", "t", "y0"],
    #                         fitalg="leastsq", weights=None)
    #fp, er, stddev = fitls(x, y, func, fp, ["a", "t", "c", "y0"],
    #                              fitalg="leastsq", weights=None)
    print fp, er
    def helper(newx):
        return func(newx, **fp)
    return helper
    