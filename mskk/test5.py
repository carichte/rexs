#!/usr/bin/env python
import os
import evaluationtools as et
import pyxrr.functions as pf
from scipy import interpolate
import pylab as pl
import mskk
import deltaf
element = "Mo"

w = 1.
dE = 1.
ns = 10
Emax = 19999.500274685415

do_interp=1

if do_interp:
    E = pl.append(Emax - pl.arange(0, 2000, dE), Emax + pl.arange(0, 2000, dE)[1:])
    print len(E)
    E.sort()
    f1, f2 = deltaf.getf(element, E, w)
    ind = pl.where(E==Emax)[0].item()
    f2[ind] = pl.mean((f2[ind-1], f2[ind+1]))
    pl.plot(E, f1, "--", E, f2, "--")
    
    E = pl.append(Emax - pl.arange(0, 2000, dE/ns), Emax + pl.arange(0, 2000, dE/ns)[1:])
    E.sort()
    f1, f2 = deltaf.getf(element, E)
    ind = pl.where(E==Emax)[0].item()
    f2[ind] = pl.mean((f2[ind-1], f2[ind+1]))
    pl.plot(E, f1, E, f2)
    f1 = et.lorentzian_filter1d(f1, w/dE*ns)
    f2 = et.lorentzian_filter1d(f2, w/dE*ns)
    #pl.plot(E, f1, E, f2)
    f1func = interpolate.interp1d(E, f1, kind="linear", bounds_error=False)
    f2func = interpolate.interp1d(E, f2, kind="linear", bounds_error=False)
    E = pl.append(Emax - pl.arange(0, 2000, dE), Emax + pl.arange(0, 2000, dE)[1:])
    E.sort()
    f2 = f2func(E)
    f1 = f1func(E)
    pl.plot(E, f1, E, f2)
else:
    E = pl.append(Emax - pl.arange(0, 2000, dE), Emax + pl.arange(0, 2000, dE)[1:])
    E.sort()
    f1, f2 = deltaf.getf(element, E, w)
    ind = pl.where(E==Emax)[0].item()
    f2[ind] = pl.mean((f2[ind-1], f2[ind+1]))
    pl.plot(E, f1, E, f2, ".-")

pl.show()

newE=E

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
pl.plot(Ere2, (f2_re_corr - f1c)-peak_re, label="unshifted")
pl.legend(loc=0)
pl.title("energy step: %.2geV, subinterpol: %ix, lorentzian broadening: %.2g"%(dE, ns, w))
pl.show()
