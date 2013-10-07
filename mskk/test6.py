#!/usr/bin/env python
import os
import evaluationtools as et
import pyxrr.functions as pf
from scipy import interpolate, special
import pylab as pl
import mskk
import deltaf
element = "Mo"

w = 0
dE = 2
if 1:
    E = pl.arange(18000., 22000., dE)
    f1, f2 = deltaf.getf(element, E)
    newE = E


peak = ((special.erf((E-20000.)/1.) + 1)/2*(f2.max() - f2.min()) + f2.min())*(E/20000.)**(-3.5/2.)
pl.plot(newE, f1, ".-")
pl.plot(newE, f2, ".-")



pl.plot(newE, peak, label="peak")
Ere, peak_re = mskk.real(newE, peak-15)
pl.plot(Ere, peak_re, label="peak trans")
Ere2, peak_im = mskk.imag(Ere, peak_re)
pl.plot(Ere2, peak_im-peak_im.min(), label="peak trans*2")

pl.show()