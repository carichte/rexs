#!/usr/bin/env python
import os
import evaluationtools as et
import pyxrr.functions as pf
from scipy import interpolate
import pylab as pl
import mskk

element = "Mo"

E, f1, f2 = pf.get_f1f2_from_db(element, table='Sasaki').T

#anchors = pl.array((12000,14000,16000,18000,22000,24000,26000,28000))
anchors = pl.array((18000.,))



f1-= pf.get_element(element)[-1]
f2func = interpolate.interp1d(E, f2)
f1func = interpolate.interp1d(E, f1)
#newE = E
#diffE =pl.diff(E)
#dE = diffE[diffE>0].mean()
dE = 1
newE = pl.arange(17000, 22000, dE)
f2 = f2func(newE)
f1 = f1func(newE)
#f1kk = mskk.real(newE, f2)


pl.plot(newE, f1, linewidth=1.5)
pl.plot(newE, f2)

Ekk, f2kk = mskk.imag(newE, f1)
pl.plot(Ekk, f2kk, label="without anchors")
Ekk, f2kk = mskk.imag(newE, f1, anc_om=anchors, anc_im=f2func(anchors))
pl.plot(Ekk, f2kk, "-", label="with %i anchors"%len(anchors))
pl.plot(anchors, f2func(anchors), "ok")
pl.legend(loc=0)

pl.show()

