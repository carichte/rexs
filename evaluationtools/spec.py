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
from collections import namedtuple

DAFS_result   = namedtuple('DAFS_result', ['integrated', 'complete', 'hkl'])
pscan_result = namedtuple('pscan_result',
                         ['time', 'voltage', 'counts', 'period', 'bins',
                          'interval', 'offset', 'duty', 'loops', 'amplitude',
                          'number'])

class SPECfile(object):
    def __init__(self, specfile):
        if os.path.isfile(str(specfile)):
            self.specfile = spec.Specfile(specfile)
        elif hasattr(specfile, "scanno"):
            self.specfile = specfile
        else:
            raise ValueError("""Input for first argument 
                                ``specfile'' not understood""")
        num = self.specfile.list().split(":")
        num = map(int, num)
        self.first = num.pop(0)
        try:
            self.last = num.pop(0)
        except IndexError:
            self.last = self.first
        self.date = self.specfile.date()
        self.epoch = self.specfile.epoch()
        self.length = self.specfile.scanno()
        self.motornames = self.specfile.allmotors()
        #self.motorpositions = namedtuple('motorpositions', self.motornames)
        self.title = self.specfile.title()
        
        self.scanno = None
        self.scan = None
        
    def __iter__(self):
        self.scanno = self.first - 1
        return self
    
    def __len__(self):
        return self.length
    
    def __getitem__(self, scanno):
        if not isinstance(scanno, int):
            raise TypeError("scan indices must be integers, not %s"\
                                              %type(scanno).__name__)
        if scanno < self.first or scanno > self.last:
            raise IndexError("scan index out of range")
        self.scan = self.specfile.select("%i"%scanno)
        return self.scan
    
    def next(self):
        if self.scanno < self.last:
            self.scanno += 1
            return self.__getitem__(self.scanno)
        else: 
            raise StopIteration


def fetch_dafs_scan(specfile, scanname):
    specfile = SPECfile(specfile)
    
    last_dafs_loop_nr = 0
    result = []
    scannum = []
    for scan in specfile:
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
            scannum.append(specfile.scanno)
    
    if len(scannum)>1:
        result.append(scannum)
    if len(result)==1:
        result = result[0]
    return result

def list_scans(specfile):
    """
        Returns a list of 2-tuples each containing the unique scan names and
        the number of scans having this name.
    """
    specfile = SPECfile(specfile)
    scannames = []
    for scan in specfile:
        comment = scan.header("Energy")
        if comment:
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


def fetch_switching_cycles(specfile):
    specfile = SPECfile(specfile)
    motors = specfile.motornames
    result = []
    scannum = []
    lastcommand = []
    lastvoltage = None
    for scan in specfile:
        command = scan.header("S")[0].split()[2:]
        positions = map(lambda s: s.partition(" ")[2], scan.header("P"))
        positions = " ".join(positions)
        positions = positions.split()
        motorpositions = dict(zip(motors, positions))
        voltage = float(motorpositions["AgilentSrc"])
        if command == lastcommand and voltage==-lastvoltage and "nix" in command:
            scannum.append(specfile.scanno)
            lastvoltage = voltage
            #if len(scannum)>1:
            #    result.append(scannum)
            #scannum = []
            continue
        else:
            if len(scannum)>1:
                result.append((abs(lastvoltage), scannum))
            scannum = [specfile.scanno]
            lastcommand = command
            lastvoltage = voltage
    return result


def process_dafs_scans(specf, indizes, trapez=True, deglitch=True, detectors=[],
                       getall=[], xiadir="", normalize=True, monitor="Io"):
    """
    Processes scans from specfile ``specf'' with scan numbers given
    in ``indizes''.
        -> normalize of all intensities to Monitor Io
        -> return as dictionary integrated values and optionally complete 
           maps
        
        Inputs:
        - trapez : bool
            integration with trapez ansatz instead of standard integration
        - deglitch : bool
            using median-filter for deglitching
        - detectors : list
            list of detectors to process
    """
    if not detectors:
        detectors = ["apd", "fluo0", "mca0", "Io"]
    sumdata = dict([(k,[]) for k in detectors])
    
    hkl = (0,0,0)
    if len(getall):
        alldata = dict([(k,[]) for k in (["q", "theta"] + getall)])
    else:
        alldata = {}
    Energy = []
    for i in indizes:
        scan = specf.select(str(i))
        try:
            hkl = scan.hkl()
        except:
            pass
            
        colname = scan.header("L")[0].split()[1:]
        #colname =  scan.alllabels()
        #print colname
        dat = scan.data()
        Energy.append(float(scan.header("Energy")[0].split()[1][:-1]))
        mon = dat[colname.index(monitor)]
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
            
            if "mca" in getall:
                alldata["mca"].append(mcadata.sum(0))
        
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
                if len(coldata) != len(q):
                    coldata = coldata[:len(q)]
            else:
                raise ValueError("Detector %s not found in Scan %i"%(col, i))
                
            if deglitch and not "apd" in col:
                coldata = ndimage.median_filter(coldata, 5)
                if col == "mca0":
                    coldata[[0, -1]] = np.median(coldata)
            if normalize and not col==monitor:
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
    res = DAFS_result(data, alldata, hkl)
    return res

def fetch_fast_pscan(specpath):
    if os.path.isfile(str(specpath)):
        sf = spec.Specfile(specpath)
    else:
        raise ValueError("Input for first argument ``specfile'' not understood")
    
    master, scannum = specpath.split("_fast_")
    scannum = int(scannum.strip(".spec"))
    master = spec.Specfile(master + ".spec")
    master = master.select(str(scannum))
    mdata = dict(zip(master.alllabels(), master.data()))
    try:
        first, last = map(int, sf.list().split(":"))
    except:
        first = int(sf.list())
        last = first
    num = range(first, last+1)
    alltimes = []
    allvolt = []
    allcnt = []
    loops = 0
    for i in num:
        scan = sf.select(str(i))
        #print mdata, i
        #mon = mdata["Mon5"][i]
        mon=0
        
        data = dict(zip(scan.alllabels(), scan.data()))
        times = data["time"]
        volt = data["fgen"]
        cnt  = data["adp"]#/mon
        #S 0 pscansettings 10 10 r 0.5 0 50 100 1
        pscansettings = filter(lambda s: "pscansettings" in s,
                               scan.header("S"))[0]
        pscansettings = pscansettings.split()
        while pscansettings.pop(0) != "pscansettings":
            pass
        period = 1./float(pscansettings[0])
        bins = int(pscansettings[6])
        interval = period/bins
        amplitude = float(pscansettings[3])/2.
        offset = float(pscansettings[4])
        duty = float(pscansettings[5])/100.
        loops += int(pscansettings[1])
        
        ind = abs(volt) > abs(amplitude)/2. # starting point
        times = times[ind]
        volt = volt[ind]
        cnt = cnt[ind]
        times -= times[0]
        times += interval/2.
        if i!=num[0]:
            times += alltimes[i-1][-1] + interval/2.
        alltimes.append(times)
        allvolt.append(volt)
        allcnt.append(cnt)
    
    times = np.hstack(alltimes)
    volt = np.hstack(allvolt)
    cnt = np.hstack(allcnt)
    result = pscan_result(times, volt, cnt, period, bins, interval, offset,
                          duty, loops, amplitude, master.number())
    return result
