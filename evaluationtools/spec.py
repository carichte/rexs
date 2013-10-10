#!/usr/bin/env python
#
#
# this is part of the 'evaluationtools' set written by Carsten Richter
# (carsten.richter@desy.de)
#
import numpy as np
import os

from scipy import ndimage, integrate
from PyMca import specfile as spec
from PyMca import EdfFile

class result(object):
    def __init__(self):
        self.integrated = None
        self.complete = None
        self.hkl = None


def fetch_dafs_scan(specfile, scanname):
    if os.path.isfile(str(specfile)):
        specfile = spec.Specfile(specfile)
    elif hasattr(specfile, "scanno"):
        pass
    else:
        raise ValueError("Input for first argument ``specfile'' not understood")
    
    num = specfile.list().split(":")
    num = map(int, num)
    
    last_dafs_loop_nr = 0
    result = []
    scannum = []
    for i in range(num[0], num[1]+1):
        scan = specfile.select(str(i))
        comment = scan.header("Energy")
        if len(comment):
            comment = comment[0].split()
        else:
            continue
        try:
            this_scanname = comment[comment.index("Scan-Name")+1]
        except IndexError:
            this_scanname = ""
        if not this_scanname == scanname:
            continue
        dafs_loop_nr = int(scan.header("DAFS")[0].split()[-1])
        
        if not dafs_loop_nr == (last_dafs_loop_nr + 1):
            #last_dafs_loop_nr = 0
            #if len(scannum)>1:
            #    result.append(scannum)
            #scannum = []
            continue
        else:
            last_dafs_loop_nr = dafs_loop_nr
            scannum.append(i)
    
    if len(scannum)>1:
        result.append(scannum)
    if len(result)==1:
        result = result[0]
    return result

def list_scans(specfile):
    if os.path.isfile(str(specfile)):
        specfile = spec.Specfile("sto_az.a")
    elif hasattr(specfile, "scanno"):
        pass
    else:
        raise ValueError("Input for first argument ``specfile'' not understood")
    
    num = specfile.list().split(":")
    num = map(int, num)
    
    scannames = []
    for i in range(num[0], num[1]+1):
        scan = specfile.select(str(i))
        comment = scan.header("Energy")
        if len(comment):
            comment = comment[0].split()
        else:
            continue
        try:
            this_scanname = comment[comment.index("Scan-Name")+1]
        except IndexError:
            this_scanname = ""
        scannames.append(this_scanname)
    
    names_unique = np.unique(scannames)
    lengths = map(scannames.count, names_unique)
    
    return zip(names_unique, lengths)

def process_dafs_scans(specf, indizes, trapez=True, deglitch = True, detectors = [], getall = [], xiadir="", normalize=True):
    """
        loading of scans with scan number in indizes of specfile specf
        norming of all intensities to Monitor Io
        - trapez : bool
            integration with trapez ansatz instead of standard integration
        - deglitch : bool
            using median-filter for deglitching
    """
    if not detectors:
        detectors = ["apd", "fluo0", "mca0", "Io"]
    sumdata = dict([(k,[]) for k in detectors])
    
    hkl = (0,0,0)
    if len(getall):
        alldata = dict([(k,[]) for k in (["q", "theta"] + getall)])
    else:
        alldata = {}
    col_mon = "Io"
    Energy = []
    for i in indizes:
        scan = specf.select(str(i))
        try:
            hkl = scan.hkl()
        except:
            pass
            
        colname = scan.header("L")[0].split()[1:]
        #colname =  scan.alllabels()
        dat = scan.data()
        Energy.append(float(scan.header("Energy")[0].split()[1][:-1]))
        mon = dat[colname.index(col_mon)]
        q = 4*np.pi*Energy[-1]/12.398 * np.sin(np.radians(dat[1]/2.))
        #q = 4*np.pi*Energy[-1]/12.398 * np.sin(np.radians(dat[0]))
        
        
        usemca = len(filter(lambda x: x not in colname, detectors))
        if usemca:
            xiaroi = scan.header('@XIAROI')
            if xiaroi:
                xiaroi = np.array(map(str.split, xiaroi))
                xiaroiname = list(xiaroi[:,1])
            
            xianame = scan.header('@XIAF')
            if xianame:
                xiapath = os.path.join(xiadir, xianame[0].split("/")[1])
                xiapath = xiapath.replace("XX", "00")
            
            if os.path.isfile(xiapath):
                edf = EdfFile.EdfFile(xiapath)
            else:
                raise ValueError("MCA File not found: %s"%xiapath)
            if edf.NumImages < 1:
                raise ValueError("MCA File empty: %s"%xiapath)
            mcadata = edf.GetData(0)
            
        
        if len(getall):
            alldata["q"].append(q)
            alldata["theta"].append(dat[1]/2.)
        for col in sumdata.keys():
            if col in colname:
                coldata = dat[colname.index(col)]
            elif xiaroi != [] and col in xiaroiname:
                ind = xiaroiname.index(col)
                roimin, roimax = xiaroi[ind, [3,4]].astype(int)
                coldata = mcadata[:, roimin:roimax].sum(1).astype(float)
            else:
                raise ValueError("Detector %s not found in Scan %i"%(col, i))
                
            if deglitch and not "apd" in col:
                coldata = ndimage.median_filter(coldata, 5)
                if col == "mca0":
                    coldata[[0, -1]] = np.median(coldata)
            if normalize and not col==col_mon:
                coldata /= mon
            if len(getall) and col in getall:
                alldata[col].append(coldata)
            dx = dat[0,1] - dat[0,0]
            if trapez:
                sumdata[col].append(integrate.trapz(coldata, x=q))
            else:
                sumdata[col].append(coldata.sum()*dx)
    data = {"Energy":Energy}
    data.update(sumdata)
    for key in data.keys():
        data[key] = np.array(data[key])
    if len(getall):
        for key in alldata.keys():
            alldata[key] = np.array(alldata[key])
    res = result()
    res.integrated = data
    res.complete = alldata
    res.hkl = hkl
    return res

