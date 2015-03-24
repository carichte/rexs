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
from rexs import tools
from rexs.xray import mskk

# Example how to use the Kramers Kronig Transformation

# correction if bad behavior on the borders, but may cause const. offset of result
corr=1


data = tools.loaddat("sto_f2.dat", skiprows=1)
#data[1]*=2
pl.plot(data[0], data[1], "-", label="$\\Delta f^{\\prime\\prime}$ measured")

KKtrans = mskk.Transform(data[0], "Ti", edges="K")

f1 = KKtrans.transform(f2=data[1], corr=corr)
pl.plot(data[0], f1, label="$\\Delta f^{\\prime}$ transformed")

f2 = KKtrans.transform(f1=f1, corr=corr)
pl.plot(data[0], f2, label="f2 transformed")

pl.plot(data[0], data[1]-f2, label="diff")

pl.xlabel("Energy (eV)")
pl.ylabel("f (electrons)")
pl.legend(loc=0, frameon=False, title="SrTiO$_3$ at Ti K-edge")
pl.grid(True)
pl.show()
