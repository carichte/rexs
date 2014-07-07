#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Mo 7. Jul 18:24:23 CEST 2014
# Computer: haso227r 
# System: Linux 3.13.0-30-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
import pylab as pl
import evaluationtools as et
import mskk

# Example how to use the Kramers Kronig Transformation

data = et.loaddat("sto_f2.dat")
pl.plot(data[0], data[1], label="f2 measured")

f1 = mskk.transform("Ti", data[0], f2=data[1], edge="K")
pl.plot(data[0], f1, label="f1 transformed")

f2 = mskk.transform("Ti", data[0], f1=f1, edge="K")
pl.plot(data[0], f2, label="f2 transformed")


pl.xlabel("Energy (eV)")
pl.ylabel("f (electrons)")
pl.legend(loc=0, frameon=False)
pl.grid(True)
pl.show()
