#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Fr 6. Jun 11:27:09 CEST 2014
# Computer: haso227r 
# System: Linux 3.2.0-63-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
#
# this is part of the 'evaluationtools' set written by Carsten Richter
#
import os
import numpy as np

class scan1d(object):
    """ 
    
        Class for one-dimensional generic scans. It consists of an arbitrary
        number of columns where one is defined as the independent variable.
        
    """
    def __init__(self, fields, data, independent=0, dtypes=None, units=None):
        """
            Initialize new scan1d class from given field names and data
        """
        
        assert hasattr(fields, "__iter__"), \
            "Need a sequence for argument #0 (fields)."
        assert all([isinstance(name, (str, unicode)) for name in fields])
        self.fields = map(str, fields)
        self.numcols = len(fields)
        shape = np.shape(data)
        assert len(shape) == 2, "2 dimensional data input expected."
        assert (shape[0] == self.numcols or shape[1] == self.numcols),\
            "Shape mismatch: lengths of data and fields disagree."
        
        if isinstance(data, np.ndarray) and shape[1] == self.numcols:
            data = data.T
        
        
        self.data = []
        
        for i in xrange(self.numcols):
            if hasattr(dtypes, "__iter__"):
                self.data.append(np.array(data[i], dtype = dtypes[i]))
            else:
                self.data.append(np.array(data[i], dtype = dtypes))
        
        if isinstance(independent, int):
            pass
        elif isinstance(independent, str):
            independent = self.fields.index(independent)
        else:
            raise ValueError("Input for `independent' not understood.")
        
        self.data.insert(0, self.data.pop(independent)) # bring x to front
        
        if units!=None:
            assert hasattr(units, "__iter__"), \
                "Need a sequence for argument units."
            assert all([isinstance(name, (str, unicode)) for name in units])
            assert len(units) == self.numcols, \
                "Same length of `units' and `fields' required."
            self.units = map(str, units)
        
        
        self.normalized = None
        
    def __len__(self):
        return self.fields.__len__()
    
    def __getitem__(self, indices):
        """
            Rewritten to handle columns names in FIOdata.colname
        """
        if isinstance(indices, str) and indices in self.fields:
            return self.data[self.fields.index(indices)]
        elif indices=="x":
            return self.data[0]
        else:
            return self.data[indices]
    
    __getattr__ = __getitem__
    
    def __dir__(self):
        return self.__dict__.keys() + dir(scan1d) + self.fields
    
    def normalize(self, col, action="divide"):
        """
            Normalizes all columns to a column specified by ``col`` and
            deletes column ``col``
        """
        if col in self.colname:
            col = self.colname.index(col)
        else:
            try: col = int(col)
            except:
                collow = map(str.lower, self.colname)
                thiscol = filter(lambda s: col.lower() in s, collow)
                col = collow.index(thiscol[0])
        for i in xrange(self.numcols):
            if i !=col and i!=0:
                if action=="divide":
                    self.data[i] /= self.data[:,col]
                elif action=="multiply":
                    self.data[i] *= self.data[:,col]
        
        self.normalized = col
    
    def to_dat(self, fname=None, comment="#", delimiter=" ", **kwargs):
        """
            Translates the scan1d data into numpy`s default columned data
            string format plus 1 header line and stores it into the file 
            ``FILENAME`` or returns it as string.
        """
        from StringIO import StringIO
        output = comment + delimiter.join(self.fields) + os.linesep
        outfile = StringIO()
        kwargs["delimiter"] = delimiter
        np.savetxt(outfile, pl.array(self.data).T, **kwargs)
        output += outfile.getvalue()
        outfile.close()
        if FILENAME == None:
            return output
        else: 
            fh = open(FILENAME, "w")
            fh.write(output)
            fh.close()

