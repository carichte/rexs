#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Mi 22. Okt 11:52:11 CEST 2014
# Computer: haso227r 
# System: Linux 3.13.0-37-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
import pylab as pl
from rexs import tools



def twogauss(x, cen1, cen2, ymax1, ymax2, w1, w2, offset):
    return tools.gaussian(x, cen1, ymax1, w1) \
         + tools.gaussian(x, cen2, ymax2, w2) \
         + offset



realparam = dict(cen1=-1, cen2=0.9, ymax1=2, ymax2=3, w1=0.5, w2=.8, offset=1)

x_m = pl.linspace(-5, 5, 501)
y_m = twogauss(x_m, **realparam) + 0.2*pl.rand(len(x_m)) # adding some noise
pl.plot(x_m, y_m, ".k", label="data", ms=2)

guess = dict((k,2*pl.rand()-1) for k in realparam) # random start values

fit = tools.fitls(x_m, y_m, twogauss, guess, fitalg="simplex")
pl.plot(x_m, fit.ysim, "-r", label="simplex")
print("Fit results simplex: %f"%fit.err)
print(fit.popt)


fit = tools.fitls(x_m, y_m, twogauss, guess, fitalg="leastsq")
pl.plot(x_m, fit.ysim, "-g", dashes=(5,2), label="leastsq")
print("Fit results leastsq: %f"%fit.err)
print(fit.popt)
pl.legend()
pl.show()