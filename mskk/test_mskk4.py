#!/usr/bin/env python
import os
import evaluationtools as et
import pyxrr.functions as pf
from scipy import interpolate, ndimage
import pylab as pl
import mskk

element = "Mo"

E, f1, f2 = pf.get_f1f2_from_db(element, table='Sasaki').T


anchors = pl.array((12000,14000,16000,18000,22000,24000,26000,28000))
#anchors = pl.array((8000.,22000))



f1-= pf.get_element(element)[-1]
f2func = interpolate.interp1d(E, f2)
f1func = interpolate.interp1d(E, f1)
#newE = E
#diffE =pl.diff(E)
#dE = diffE[diffE>0].mean()
dE = 5
newE = pl.arange(10000, 30000, dE)
f2 = f2func(newE)
f1 = f1func(newE)
#f1kk = mskk.real(newE, f2)


if 1:
    f1 = et.lorentzian_filter1d(f1, 10)
    f2 = et.lorentzian_filter1d(f2, 10)


pl.plot(newE, f1, linewidth=1.5)
pl.plot(newE, f2)

Ekk, f1kk = mskk.real(newE, f2)
pl.plot(Ekk, f1kk, label="without anchors")
Ekk, f1kk = mskk.real(newE, f2, anc_om=anchors, anc_re=f1func(anchors))
pl.plot(Ekk, f1kk, "-", label="with %i anchors"%len(anchors))
pl.plot(anchors, f1func(anchors), "ok")
pl.legend(loc=0)

pl.show()


pl.plot(newE, f2)
pl.plot(Ekk, f1kk, "-", label="with %i anchors"%len(anchors))

Ekk2, f2kk = mskk.imag(Ekk, f1kk)
pl.plot(Ekk2, f2kk, label="without anchors")
Ekk2, f2kk = mskk.imag(Ekk, f1kk, anc_om=anchors, anc_im=f2func(anchors))
pl.plot(Ekk2, f2kk, "-", label="with %i anchors"%len(anchors))
pl.plot(anchors, f2func(anchors), "ok")
pl.legend(loc=0)

pl.show()


f1 = ndimage.gaussian_filter1d(f1func(newE), 1)


pl.plot(newE, f2)
pl.plot(newE, f1, "-")
pl.plot(Ekk, f1kk, "-")

Ekk2, f2kk = mskk.imag(Ekk, f1kk)
pl.plot(Ekk2, f2kk, label="without anchors")
pl.legend(loc=0)

pl.show()
