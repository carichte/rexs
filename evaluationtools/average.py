#!/usr/bin/env python
#
#
# this is part of the 'evaluationtools' set written by Carsten Richter
# (carsten.richter@desy.de)
#


import numpy as np
import os
from scipy import ndimage, optimize, special
import time
from diverse import loaddat
from FIOdata import FIOdata

def reduce_measured_data(data, new_stepping):
    if len(data.shape)==1: data=np.vstack((np.arange(len(data)), data)).T
    elif data.shape[1]==1: data=np.vstack((np.arange(len(data)), data)).T
    theta_curr=data[:,0].min() + new_stepping/2. #startpunkt
    columns=""
    results=None
    column_amount=len(data[0])
    zeile=np.zeros(column_amount)
    anzahl=0
    for j in range(len(data[:,0])):
        theta=data[j,0]
        if theta<(theta_curr+new_stepping/2):
            zeile[1]+=data[j,1]
            anzahl+=1
            if column_amount>2: 
                if j>0: zeile[2]=1./np.sqrt(data[j,2]**(-2) + zeile[2]**(-2)) # 1/sigma
                else: zeile[2] = data[j,2]
        else:
            zeile[0]=theta_curr
            zeile[1]/=anzahl
            if results==None: results=zeile
            else: results=np.vstack((results, zeile))
            theta_curr+=new_stepping
            anzahl=1
            zeile=data[j].copy()
    return (results)


def average_fio_repeats(FILENAME):
    """
        Averages the y-data of a repeated measurement.
        
        Input one of the .fio filenames.
    """
    root = FILENAME.strip(".fio").rpartition("r")
    CollData = FIOdata(FILENAME)
    if root[2].isalnum: root = root[0]
    else: return CollData
    CollData.data[:,1:]*=0.
    FILEDIR=os.path.dirname(FILENAME)
    if not FILEDIR: FILEDIR = "."
    repeat = 0
    for fname in os.listdir(FILEDIR):
        if root not in fname or ".fio" not in fname: continue
        CollData.data[:,1:] += FIOdata(os.path.join(FILEDIR, fname)).data[:,1:]
        repeat+=1
    CollData.repeats = repeat
    CollData.comment += "Repeats: %i\n" %repeat
    return CollData
    

def average_fio_files(include, try_repair=False):
    """
        Averages the y-data of the set of .fio files containing the string *include*.
        
        In Continous scans there is sometimes a length mismatch of 1.
        Therefore it is tried to repair this by appending the last date of the 
        shorter scan, if try_repair == True.
    """
    CollData = None
    FILEDIR=os.path.dirname(include)
    if not FILEDIR: FILEDIR = "."
    for fioname in os.listdir(FILEDIR):
        if include not in fioname or ".fio" not in fioname: continue
        if not CollData: CollData = FIOdata(FILEDIR + "/" + fioname)
        else:
            NextData = FIOdata(FILEDIR + "/" + fioname)
            if len(NextData.data) == len(CollData.data): pass
            elif not try_repair: raise ValueError("Measurements differ in number of data-points")
            elif len(NextData.data) == len(CollData.data)+1: CollData.data = np.vstack((CollData.data, CollData.data[-1]))
            elif len(NextData.data) == len(CollData.data)-1:  NextData.data = np.vstack((NextData.data, NextData.data[-1]))
            else: raise ValueError("Measurements differ in number of data-points")
            CollData.data[:,1:] += NextData.data[:,1:]
    return CollData
    


def average_dat(include, exclude=[], path=os.curdir, yerr=False, sumonly=False, col_weights=None):
    """
        Average columned data files in a certain path.
        Including and excluding certain substrings of the filename.
        
        1st row can be comments.
        1st column = x
        assumes data to be standard deviations if yerr==True
    """
    flist = sorted(os.listdir(path))
    #flist.reverse()
    num = None
    for fname in flist:
        skipfile = False
        if hasattr(include, "__iter__"):
            for tag in include:
                skipfile += (tag not in fname)
        elif isinstance(include, str):
            if include not in fname: continue
        else: raise ValueError("Data Type of include argument not understood")
        if hasattr(exclude, "__iter__"):
            for tag in exclude: 
                skipfile += (tag in fname)
        elif isinstance(exclude, str):
            if exclude in fname: continue
        else: raise ValueError("Data Type of exclude argument not understood")
        if skipfile: continue
        
        fpath = os.path.join(path, fname)
        print("Opening %s..." %fname)
        if num==None:
            data, header = loaddat(fpath)
            num = np.ones(len(data))
            if col_weights!=None:
                try: col_weights = int(col_weights)
                except: col_weights = header.split().index(col_weights)
                for col in np.arange(1, len(data[0])):
                    if col!=col_weights:
                        data[:,col] *= data[:,col_weights]
        else:
            try: nextdata = np.loadtxt(fpath, skiprows=0)
            except: nextdata = np.loadtxt(fpath, skiprows=1)
            if col_weights!=None:
                for col in np.arange(1, len(data[0])):
                    if col!=col_weights:
                        nextdata[:,col] *= nextdata[:,col_weights]
            ld = len(data)
            ln = len(nextdata)
            ind = min(ld,ln)
            if yerr:
                data[:ind,1:] = np.sqrt(data[:ind,1:]**2 + nextdata[:ind,1:]**2)
            else:
                data[:ind,1:]+=nextdata[:ind,1:]
            num[:ind]+=1
            if ln>ld:
                data = np.vstack((data, nextdata[ind:]))
                num = np.append(num, np.ones(ln-ld))
    if num==None: return 0
    if not sumonly:
        if col_weights == None:
            data[:,1:]/=num[:,np.newaxis]
        else:
            for col in np.arange(1,len(data[0])):
                if col!=col_weights:
                    data[:,col]/=data[:,col_weights]
    return data, header, num

def average_dat_xye(include, exclude=chr(0), path=os.curdir):
    """
        Average columned data files in a certain path.
        Including and excluding certain substrings of the filename.
        
        1st row can be comments.
        1st column = x
        even columns = y
        odd columns = y_err
        (x y yerr y yerr ... )
    """
    flist = os.listdir(path)
    num = 0
    for fname in flist:
        if include not in fname or exclude in fname: continue
        print fname
        fpath = os.path.join(path, fname)
        if not num:
            data, header = loaddat(fpath)
            ind = ( np.arange(len(data[0]))%2 ) == 1
            errind = ( np.arange(len(data[0]))%2 ) == 0
            errind[0]=False
        else:
            try: nextdata = np.loadtxt(fpath, skiprows=0)
            except: nextdata = np.loadtxt(fpath, skiprows=1)
            data[:,ind]+=nextdata[:,ind]
            data[:,errind] = np.sqrt(data[:,errind]**2 + nextdata[:,errind]**2)
        num+=1
    #if not sumonly: data[:,1:]/=num
    return data, header


def rebin(x, y, weights=None, bins=None, xmin=None, xmax=None,
          discard_empty=False):
    """
    
        Function that averages subsequent datasets via histogramming. The data
        does not have to be present for repeatingly the same values of the
        independent variable and does not have to be ordered. Therefore, a
        rebinning to a new equi-distant x-axis takes place. 
        
    """
    x = np.hstack(x)
    if xmin==None:
        xmin = x.min()
    if xmax==None:
        xmax = x.max()
    ind = (x>=xmin) * (x<=xmax)
    x = x[ind]
    y = np.hstack(y)[ind]
    if bins==None:
        bins = (x.max()-x.min())/np.diff(np.sort(x)).max()
        bins = int(np.floor(bins))
    if weights==None:
        weights = np.ones(len(x))
    else:
        weights = np.hstack(weights)[ind]
    y, newx = np.histogram(x, bins=bins, weights=y)
    num, newx = np.histogram(x, bins=bins, weights=weights)
    dx = newx[1] - newx[0]
    x = newx[1:] - dx/2.
    y /= num
    ind = num > 0
    if not ind.all() and discard_empty:
        y = y[ind]
        x = x[ind]
    return x, y

