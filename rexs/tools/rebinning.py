import numpy as np

def rebin(x, y, weights=None, bins=None, xmin=None, xmax=None,
          discard_empty=False):
    """
    
        Function that averages subsequent datasets via histogramming. The data
        does not have to be present for repeatingly the same values of the
        independent variable and does not have to be ordered. Therefore, a
        rebinning to a new equi-distant x-axis takes place. 
        
    """
    #x = np.hstack(x)
    if np.ndim(x) > 1:
        x = np.ravel(x)
        y = np.ravel(y)
    if xmin is None:
        xmin = x.min()
    if xmax is None:
        xmax = x.max()
    ind = (x>=xmin) * (x<=xmax)
    x = x[ind]
    y = y[ind]
    if bins is None:
        bins = (x.max()-x.min())/np.diff(np.sort(x)).max()
        bins = int(np.floor(bins))
    if weights is None:
        weights = np.ones(len(x))
    else:
        weights = np.ravel(weights)[ind]
    y, newx = np.histogram(x, bins=bins, weights=y, range=(xmin, xmax))
    num, newx = np.histogram(x, bins=bins, weights=weights, range=(xmin, xmax))
    dx = newx[1] - newx[0]
    x = newx[1:] - dx/2.
    y /= num
    ind = num > 0
    if not ind.all() and discard_empty:
        y = y[ind]
        x = x[ind]
    return x, y



