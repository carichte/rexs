#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Do 18. Dez 18:42:23 CET 2014
# Computer: haso227r 
# System: Linux 3.13.0-40-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter All rights reserved.
#----------------------------------------------------------------------
"""
    Simple text-based selection of scans and columns from `.spec` file for
    plotting using `matplotlib`.
"""
import os
import sys 
import argparse
import pylab as pl
import evaluationtools
from evaluationtools._custom_screen import Screen
from evaluationtools import spec
import signal
import ConfigParser



parser = argparse.ArgumentParser(description=
   "Simple text-based selection of scans and columns from `.spec` file for "
   "plotting using `matplotlib`.")
#parser.add_argument("-s", "--sort", choices=["time", "name"], default=None,
#                    help="sort input files")
parser.add_argument("specfile", nargs=1, default=[], 
                    help=".spec file to open")
parser.add_argument("scanno", type=int, default=None,
                    help="Numbers of scans to plot", nargs="*")
parser.add_argument("-c", "--col", type=str, default=None,
                    help="Name of columns to plot for each scan")
#parser.add_argument("-o", "--outfile", type=str, 
#                    default="fdmnes_overview.csv", help="output file")
args = parser.parse_args()

fname = args.specfile[0]

if not os.path.isfile(fname):
    raise ValueError("File not found: %s"%fname)

sf = spec.SPECfile(fname)




confdir = os.path.join(os.path.expanduser("~"), ".specviewer")
if not os.path.isdir(confdir):
    try:
        os.mkdir(confdir)
    except:
        pass
conffile = os.path.join(confdir, "cache")
if not os.path.isfile(conffile):
    try:
        open(conffile, 'w').close()
    except:
        pass

# READING CONFIG ##########################################################
if os.path.isfile(conffile):
    conf = ConfigParser.ConfigParser()
    conf.read(conffile)
    if not conf.has_option(fname, "plotscans"):
        plotscans = []
    else:
        plotscans = conf.get(fname, "plotscans").split()
    
    if not conf.has_option(fname, "cols"):
        cols = []
    else:
        cols = conf.get(fname, "cols").split()
else:
    plotscans, cols = [], []

# CHOOSING SCANS TO PLOT ##################################################
if not args.scanno:
    scaninfo = [(str(scan.number()), scan.command()) for scan in sf]
    numbers, commands = zip(*scaninfo)
    scaninfo = dict(scaninfo)
    with Screen() as screen:
        signal.signal(signal.SIGWINCH, screen.resize)
        screen.clear()
        screen.draw_title("Use: <space> - mark scan to plot |"
                          " <enter> - plot scan | ENTER - plot/proceed")
        
        pos = plotscans[0] if plotscans else 0
        attr = screen.menu(map(str,numbers), pos, info=scaninfo, showall=True,
                           selected=plotscans, exitbutton=ord("\n"))
        
        if attr=="exit":
            pass
        else:
            plotscans = [int(attr)]
    plotscans = map(int, plotscans)
else:
    plotscans = args.scanno

if not plotscans:
    raise StandardError("No scan chosen.")

# CHOOSING COLUMNS TO PLOT ################################################
if args.col == None:
    labels = set(sum([sf[i].alllabels()[1:] for i in plotscans], []))
    with Screen() as screen:
        signal.signal(signal.SIGWINCH, screen.resize)
        screen.clear()
        screen.draw_title("Use: <space> - mark columns to plot |"
                          " <enter> - plot colum | p - plot/proceed")
        
        pos = cols[0] if cols else 0
        attr = screen.menu(sorted(labels), pos, showall=True, 
                           selected=cols, debug=True)
        
        if attr=="exit":
            pass
        else:
            cols = [attr]
else:
    cols = [args.col]

if not cols:
    raise StandardError("No columns chosen.")



# PLOTTING ################################################################
nrow = int(pl.sqrt(len(cols)))
ncol = int(pl.ceil(float(len(cols))/nrow))
fig, ax = pl.subplots(nrow, ncol, sharex=True, squeeze=False)

xnames = set()
for i in plotscans:
    scan = sf[i]
    labels = scan.alllabels()
    xnames.add(labels[0])
    for j, colname in enumerate(cols):
        if not colname in labels:
            raise ValueError("Column `%s` not found."%colname)
        ax.ravel()[j].plot(scan.datacol(1), scan.datacol(colname), 
                           label="#%i:  %s"%(i, scan.command()))

for j, colname in enumerate(cols):
    ax.ravel()[j].set_xlabel(r" / ".join(xnames))
    ax.ravel()[j].set_ylabel(colname)
pl.legend()
pl.show()


# STORING CURRENT SELECTION ###############################################
if os.path.isfile(conffile):
    conf.read(conffile)
    if not conf.has_section(fname):
        conf.add_section(fname)
    conf.set(fname, "plotscans", " ".join(map(str, plotscans)))
    conf.set(fname, "cols", " ".join(cols))
    with open(conffile, "w") as fh:
        conf.write(fh)

