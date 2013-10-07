#!/usr/bin/env python
import os
import evaluationtools as et
import pyxrr.functions as pf
from scipy import interpolate
import pylab as pl
import mskk

element = "Nb"

E, f1, f2 = pf.get_f1f2_from_db(element, table='Sasaki').T



f1-= pf.get_element(element)[-1]
f2func = interpolate.interp1d(E, f2)
f1func = interpolate.interp1d(E, f1)

dE = 5
newE = pl.arange(15000, 25000, dE)
f2 = f2func(newE)
f1 = f1func(newE)


pl.plot(newE, f2)
pl.plot(newE, f1, "-")


Ekk2, f2kk = mskk.imag(newE, f1)
pl.plot(Ekk2, f2kk, label="without anchors")
pl.legend(loc=0)

pl.show()
