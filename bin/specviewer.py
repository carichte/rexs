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
import hashlib
from evaluationtools._custom_screen import Screen
from evaluationtools import spec
import evaluationtools as et
import signal
import ConfigParser



#print __file__

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

fname = os.path.abspath(args.specfile[0])

if not os.path.isfile(fname):
    raise ValueError("File not found: %s"%fname)

sf = spec.SPECfile(fname)




confdir = os.path.join(os.path.expanduser("~"), ".cache", "specviewer")
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

def md5_for_file(path, block_size=256*128, hr=False):
    '''
    Block size directly depends on the block size of your filesystem
    to avoid performances issues
    Here I have blocks of 4096 octets (Default NTFS)
    '''
    md5 = hashlib.md5()
    with open(path,'rb') as f: 
        for chunk in iter(lambda: f.read(block_size), b''): 
             md5.update(chunk)
    if hr:
        return md5.hexdigest()
    return md5.digest()

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
        submitoptions = {"proceed":ord("p")}
        attr = screen.menu(map(str,numbers), pos, submit=submitoptions, 
                           info=scaninfo, showall=True, selected=plotscans)
        
        if attr in ["proceed", "__defaultaction"]:
            pass
        elif attr=="cancel":
            raise RuntimeError("Aborted by user.")
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
    nbmca = min([sf[i].nbmca()/sf[i].lines() for i in plotscans])
    for i in xrange(nbmca):
        labels.add("__MCA_%i_3d"%(i+1))
        labels.add("__MCA_%i_3d-log"%(i+1))
        labels.add("__MCA_%i_2d-sum"%(i+1))
        labels.add("__MCA_%i_2d-roi"%(i+1))
    with Screen() as screen:
        signal.signal(signal.SIGWINCH, screen.resize)
        screen.clear()
        screen.draw_title("Use: <space> - mark columns to plot |"
                          " <enter> - plot colum | p - plot/proceed")
        
        pos = cols[0] if cols else 0
        labels = sorted(labels)
        submitoptions = {"save to .dat":ord("s"), "plot":ord("p")}
        attr = screen.menu(labels, pos, showall=True, submit=submitoptions,
                           selected=cols, debug=True)
        
        #raise
        saveonly = False
        if attr in ["plot", "__defaultaction"]:
            pass
        elif attr=="save to .dat":
            saveonly = True
        elif attr=="cancel":
            raise RuntimeError("Aborted by user.")
        else:
            cols = [attr]
else:
    cols = [args.col]

if not cols:
    raise StandardError("No columns chosen.")


print(attr)

# PLOTTING ################################################################
nrow = int(pl.sqrt(len(cols)))
ncol = int(pl.ceil(float(len(cols))/nrow))
if not saveonly:
    fig, ax = pl.subplots(nrow, ncol, sharex=True, squeeze=False)
    xnames = set()
    vmin = ncol * nrow * [pl.inf]
    vmax = ncol * nrow * [-pl.inf]
else:
    prefix = et.ask("Enter prefix for data file output", "Scan")


if filter(lambda s: "__MCA_" in s and "roi" in s, cols):
    nint = et.ask("Number of channels to integrate: +-", 10, int)
    doBGcorrect = et.yesno("Perform 2nd order polynomial background subtraction?")

for i in plotscans:
    scan = sf[i]
    x = scan.datacol(1)
    labels = scan.alllabels()
    if saveonly:
        fout = r"%s_%i.dat"%(prefix, i)
        savepath = os.path.join(os.path.splitext(fname)[0], fout)
    
        savecols = [x]
        header = [labels[0]]
    else:
        xnames.add(labels[0])
    
    print("Plotting scan: %s ..."%("#%i:  %s"%(i, scan.command())))
    for j, colname in enumerate(cols):
        if not saveonly:
            thisax = ax.ravel()[j]
        if colname in labels:
            if not saveonly:
                thisax.plot(x, scan.datacol(colname), 
                                   label="#%i:  %s"%(i, scan.command()))
                lblax = thisax
            else:
                savecols.append(scan.datacol(colname))
                header.append(colname)
        elif colname.startswith("__MCA_"):
            i_mca = int(colname.rsplit("_",2)[1])
            anf = scan.lines() * (i_mca-1)
            end = scan.lines() * i_mca
            data = [scan.mca(line+1) for line in xrange(anf, end)]
            data = pl.vstack(data).T
            channels = pl.arange(data.shape[0])
            if doBGcorrect:
                for l in xrange(data.shape[1]):
                    indf = data[:,l] < (pl.median(data[:,l]))
                    poly = et.PolynomialFit(channels, data[:,l], indf=indf)
                    data[:,l] -= poly
            if "3d" in colname and not saveonly:
                if "-log" in colname:
                    imdata = pl.log(data)
                    ind = pl.isfinite(imdata)
                else:
                    imdata = data
                    ind = slice(None)
                vmin[j] = min(vmin[j], imdata[ind].min())
                vmax[j] = max(vmax[j], imdata[ind].max())
                thisax.imshow(pl.flipud(imdata), aspect="auto",
                                 vmin=vmin[j], vmax=vmax[j],
                                 extent=(x[0], x[-1], 0, data.shape[0]-1))
                thisax.set_title("#%i:  %s"%(i, scan.command()))
                thisax.grid(False)
            elif "2d" in colname:
                if "roi" in colname:
                    imax = data.sum(1).argmax()
                    ind = slice(imax - nint, imax + nint)
                    print("Channels in roi: %i...%i"%(imax - nint, imax + nint))
                elif "sum" in colname:
                    ind = slice(None)
                if not saveonly:
                    thisax.plot(x, data[ind].sum(0), 
                                   label="#%i:  %s"%(i, scan.command()))
                    lblax = thisax
                else:
                    savecols.append(data[ind].sum(0))
                    header.append(colname)
                
                
            
        else:
            raise ValueError("Column `%s` not found."%colname)
    
    if saveonly:
        if not os.path.isdir(os.path.dirname(savepath)):
            os.mkdir(os.path.dirname(savepath))
        et.savedat(savepath, savecols, " ".join(header))

if not saveonly:
    for j, colname in enumerate(cols):
        ax.ravel()[j].set_xlabel(r" / ".join(xnames))
        ax.ravel()[j].set_ylabel(colname)
    lblax.legend()
    pl.show()

print vmin,vmax
# STORING CURRENT SELECTION ###############################################
if os.path.isfile(conffile):
    conf.read(conffile)
    if not conf.has_section(fname):
        conf.add_section(fname)
    conf.set(fname, "plotscans", " ".join(map(str, plotscans)))
    conf.set(fname, "cols", " ".join(cols))
    with open(conffile, "w") as fh:
        conf.write(fh)

