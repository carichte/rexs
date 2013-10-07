#!/usr/bin/env python
import os
import evaluationtools as et
import pyxrr.functions as pf
from scipy import interpolate
import pylab as pl
import mskk
import deltaf
element = "Mo"

w = 1


dE = .5
do_interp=1
#lw=60.
Erange = 2000.
lw=min(Erange/w*2, 200.)

if do_interp:
    E, ind = deltaf.get_energies(element, 17000, 23000, fwhm_ev=w)
    print len(E)
    f1 = deltaf.getfquad(element, E, w, f1f2='f1', lw=lw)
    f2 = deltaf.getfquad(element, E, w, f1f2='f2', lw=lw)
    
    if w==0:
        f2[ind] = pl.mean((f2[ind-1], f2[ind+1]))
    pl.plot(E, f1, "--", E, f2, "--")
    
    f1func = interpolate.UnivariateSpline(E, f1, k=3, s=0)
    f2func = interpolate.UnivariateSpline(E, f2, k=3, s=0)
    #f2func = interpolate.interp1d(E, f2, kind="linear", bounds_error=False)
    E = pl.append(E[ind] - pl.arange(0, Erange, dE), E[ind] + pl.arange(0, Erange, dE)[1:])
    print len(E)
    E.sort()
    f2 = f2func(E)
    f1 = f1func(E)
    pl.plot(E, f1, E, f2)
else:
    E = pl.append(Emax - pl.arange(0, Erange, dE), Emax + pl.arange(0, Erange, dE)[1:])
    E.sort()
    f1, f2 = deltaf.getf(element, E, w)
    ind = pl.where(E==Emax)[0].item()
    f2[ind] = pl.mean((f2[ind-1], f2[ind+1]))
    pl.plot(E, f1, E, f2, ".-")

pl.show()

newE = E

peak = et.gaussian(newE, 20000., 5., 100, 0)

pl.plot(newE, f2, label="f2")
f2 += peak

pl.plot(newE, f2, label="f2+peak")
pl.plot(newE, peak, label="peak")

Ere, peak_re = mskk.real(newE, peak)
Ere2, f2_re = mskk.real(newE, f2, corr=False)
Ere2, f2_re_corr = mskk.real(newE, f2)


pl.plot(Ere, peak_re, label="peak trans")
pl.plot(Ere2, f2_re_corr, label="f2 trans corr")

if do_interp:
    f1c = f1func(Ere2)
else:
    f1c, f2c = deltaf.getf(element, Ere2, w)
pl.plot(Ere2, f1c, ".-", label="f1 tab shift")
pl.plot(Ere2, f2_re - f1c, label="f2 trans - f1")
pl.plot(Ere2, f2_re_corr - f1c, label="f2 trans corr - f1")


pl.legend()

pl.figure()
#pl.subplot(121)
#pl.plot(Ere2, (f2_re_corr - f1c)-peak_re, label="unshifted")
pl.subplot(111)
pl.plot(Ere2[1:], pl.diff((f2_re_corr - f1c)-peak_re), label="unshifted")
pl.legend(loc=0)
pl.ylim(ymin=-1e-4, ymax=1e-4)
pl.title("energy step: %.2geV, lorentzian broadening: %.2g, lorentz integration width: %g"%(dE, w, lw))
pl.show()
