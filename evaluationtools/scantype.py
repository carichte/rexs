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
import diverse
import re

class scan1d(object):
    """ 
    
        Class for one-dimensional generic scans. It consists of an arbitrary
        number of columns where one is defined as the independent variable.
        
    """
    def __init__(self, data, fields=None, independent=0, dtypes=None, units=None):
        """
            Initialize new scan1d class from given field names and data
        """
        
        if fields==None:
            if hasattr(data, "alllabels") and hasattr(data, "data"):
                fields = data.alllabels()
                data = data.data()
            elif isinstance(data, (str, unicode)) and os.path.isfile(data):
                data, fields = diverse.loaddat(data, parse=False)
                fields = fields.split()
            else:
                raise ValueError("Need input for 2nd argument (`fields').")
        assert hasattr(fields, "__iter__"), \
            "Need a sequence for argument #0 (fields)."
        assert all([isinstance(name, (str, unicode)) for name in fields])
        units = []
        for i, field in enumerate(fields):
            field = str(field)
            unit  = re.match("_\((.+)\)",field)
            if unit != None:
                units.append(unit.groups()[-1])
            else:
                units.append("")
            
            field = re.match("([a-zA-Z][a-zA-Z0-9]*)",field)
            if field==None:
                fields[i] = "field%i"%i
            else:
                fields[i] = field.group(0)
        self.fields = fields
        self.units = units
        
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
        
    def from_specscan(scan):
        return scan1d(scan.alllabels(), scan.data())
    
    def __len__(self):
        return self.numcols
    
    def __setitem__(self, attr, val):
        """
            Rewritten to handle columns names in scan1d.fields
        """
        if isinstance(attr, str) and attr in self.fields:
            self.data[self.fields.index(attr)] = val
        elif attr=="x":
            self.data[0] = val
        elif isinstance(attr, (int, long)):
            self.data[attr] = val
        else:
            self.__dict__[attr] = val
    
    def __getitem__(self, attr):
        """
            Rewritten to handle columns names in scan1d.fields.
            Only called when `attr' is actually not found.
        """
        if attr=="fields":
            return []
        elif isinstance(attr, str) and attr in self.fields:
            return self.data[self.fields.index(attr)]
        elif attr=="x":
            return self.data[0]
        elif isinstance(attr, int):
            return self.data[attr]
        else:
            raise ValueError("Element not found: %s"%str(attr))
    

    __getattr__ = __getitem__
    __setattr__ = __setitem__
    
    def __dir__(self):
        return self.__dict__.keys() + dir(scan1d) + self.fields
    
    def crop(self, ind):
        """
            Crops points of the scan corresponding to the 1d-array
            `ind' which contains a boolean array mask.
            
            Input:
                ind : numpy.ndarray, dtype=bool, ndim=1
                    defines which points to take.
                    The length has to agree with the length of the scan.
        """
        assert (ind.size==self.data[0].size), "Invalid length of `ind'."
        for i in xrange(self.numcols):
            self.data[i] = self.data[i][ind]
    
    def normalize(self, col, action="divide"):
        """
            Normalizes all columns to a column specified by ``col`` and
            deletes column ``col``
        """
        if col in self.fields:
            col = self.fields.index(col)
        else:
            try: col = int(col)
            except:
                collow = map(str.lower, self.fields)
                thiscol = filter(lambda s: col.lower() in s, collow)
                col = collow.index(thiscol[0])
        for i in xrange(self.numcols):
            if i !=col and i!=0:
                if action=="divide":
                    self.data[i] /= self.data[col]
                elif action=="multiply":
                    self.data[i] *= self.data[col]
        
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
        np.savetxt(outfile, np.array(self.data).T, **kwargs)
        output += outfile.getvalue()
        outfile.close()
        if fname == None:
            return output
        else: 
            fh = open(fname, "w")
            fh.write(output)
            fh.close()

