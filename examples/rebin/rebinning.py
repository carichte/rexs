#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Mi 17. Sep 17:40:07 CEST 2014
# Computer: haso227r 
# System: Linux 3.13.0-35-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------

if __name__=="__main__":
    import pylab as pl
    from rebin import rebin1d
    x = pl.linspace(-20,20,1001)
    num = pl.ones(1001)
    y = pl.sin(x) + 0.*pl.rand(len(x))
    newx = pl.linspace(-20,20,101)
    newy = rebin1d(x,y,newx, average=True)
    
    pl.plot(x,y, "-", label="original")
    pl.plot(newx,newy, ".-", label="averaged")
    pl.legend(frameon=True)
    
    
    
    a = pl.rand(1001)
    y, bins = pl.histogram(a, bins=100)
    x = bins[1:] - pl.diff(bins)[0]/2.
    pl.figure()
    pl.plot(x,y, label="original")
    newbins = pl.linspace(0,1,21)
    newy = rebin1d(bins, y, newbins, average=True)
    newx = newbins[1:] - pl.diff(newbins)[0]/2.
    pl.plot(newx,newy, "o-", label="averaged")
    newy = rebin1d(x, y, newx, average=True)
    pl.plot(newx,newy, "x-", label="averaged 2")
    pl.legend(loc=0, frameon=True)
    pl.show()


