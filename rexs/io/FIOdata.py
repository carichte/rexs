#
#
# this is part of the rexs package written by Carsten Richter
# (carsten.richter@desy.de)
#

from __future__ import print_function
import os
import sys
import time
import locale

if sys.version_info[0]<3:
    from StringIO import StringIO
else:
    from io import StringIO

import numpy as np

class FIOdata(object):
    """ 
        This class handles measurement data files that are present in
        the .fio format which is produced at the DESY Photon Science
        Instruments.
    """
    def __init__(self, FILENAME, verbose=False):
        """
            This opens a .fio file using a path to the file or a file
            handle. If verbose is True, additional information is printed.
            
            After initialization the following objects are available:
            .name       - name of measurement
            .repeats    - number of repeat
            .parameters - dictionary of motor positions
            .colname    - list of columns names
            .comment    - comment string
            .data       - numpy array of measured data
            .startsec   - unix time of measurement start
            .stopsec    - unix time of measurement stop
            .starttime  - struct_time of measurement start
            .stoptime   - struct_time of measurement stop
            
            Methods:
            .normalize  - method to normalize all columns onto a selected one
            .to_dat     - convert to .dat format (columns with 1 header line)
        """
        if isinstance(FILENAME, str):
            if verbose: print("Loading %s"%FILENAME)
            data=open(FILENAME, "r")
        elif hasattr(FILENAME, "readline"):
            if verbose: print("Loading %s"%FILENAME.name)
            data = FILENAME
        else: raise ValueError('fname must be a string or file handle')
        
        self.comment=""
        self.parameters={}
        colname=[]
        self.repeats = 1
        #self.data=np.array([])
        flag = True
        while flag:
            line = data.readline()
            if not line: break
            if not line.find("!"):
                continue
            if "%c" in line:
                line = data.readline()
                while not line.startswith("!"):
                    if not line: break
                    self.comment+=line
                    line = data.readline()
            if "%p" in line:
                line = data.readline()
                while not ("! Data" in line):
                    if line.startswith("!"): 
                        line = data.readline()
                        continue
                    elif not line: break
                    line=line.replace(" ","")
                    [param, value]=line.split("=")
                    self.parameters[param]=float(value)
                    line = data.readline()
            if "%d" in line:
                line = data.readline()
                while not line.startswith("!"):
                    if not line: break
                    if "Col" in line:
                        colname.append(line.split()[2])
                    else: 
                        flag = False
                        break
                    line = data.readline()
        numdata = StringIO(line + data.read())
        data.close()
        self.data=np.loadtxt(numdata, comments="!")
        i=0
        cond = True
        if len(colname)<=1:
            self.name = colname[0]
            self.colname = colname
        else:
            while cond:
                i+=1
                for name in colname:
                    cond *= (name[:i]==colname[0][:i])
            self.name = colname[0][:(i-2)]
            for j in range(len(colname)):
                colname[j] = colname[j][(i-1):]
            self.colname = colname
        
        words = self.comment.split()
        
        if "sampling" in words:
            ind = words.index("sampling")
            self.sampletime = float(words[ind+1])
        try:
            i1 = words.index("ended")
            day = words[i1-2]
            time0 = words[i1-1][:-1]
            timeE = words[i1+1]
            orglocal = locale.getlocale(locale.LC_TIME)
            locale.setlocale(locale.LC_TIME, ("en","UTF8"))
            self.starttime = time.strptime(day + time0, "%d-%b-%Y%H:%M:%S")
            self.stoptime = time.strptime(day + timeE, "%d-%b-%Y%H:%M:%S")
            locale.setlocale(locale.LC_TIME, orglocal)
            self.startsec = time.mktime(self.starttime)
            self.stopsec = time.mktime(self.stoptime)
        except Exception as err:
            print("Warning: starting or stopping time of scan could not be",
                  "determined: ", err)
            self.starttime = np.nan
            self.stoptime = np.nan
            self.startsec = np.nan
            self.stopsec = np.nan
        
    def __len__(self):
        return self.data.__len__()
    
    def __getitem__(self, indices):
        """
            Rewritten to handle columns names in FIOdata.colname
        """
        if isinstance(indices, str) and indices in self.colname:
            return self.data[:,self.colname.index(indices)]
        else:
            return self.data[indices]
    def __repr__(self):
        return self.comment
    def parameters_nice(self, format="%.4g"):
        """
            This function just returns a string containing the motor position
            and further parameters as they were during the measurement.
            It is preformatted in a well readable way.
            The ``format`` of the floating point numbers  can be specified
            using the common string formatting presentation types.
            
            The bare information is located in the FIOdata.parameters
            dictionary.
        """
        s=""
        p = self.parameters
        length = max([len(key) + len(format%val) 
                      for (key,val) in p.iteritems()])
        length += 2
        for key in sorted(p):
            s += key + (length-len(key)-len(format%(p[key]))) * " " \
                     +  format%(p[key]) + os.linesep
        
        
        
        return s
    def normalize(self, col):
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
        for i in range(len(self.colname)):
            if i !=col and i!=0:
                self.data[:,i] /= self.data[:,col]
        ind = np.ones(len(self.colname)).astype(bool)
        ind[col]=False
        self.data = self.data[:,ind]
        self.colname.pop(col)
    def to_dat(self, FILENAME=None, delimiter=" ", fmt="%.18e"):
        """
            Translates the FIOdata.data into numpy`s default columned data
            string format plus 1 header line and stores it into the file 
            ``FILENAME`` or returns it as string.
        """
        output = delimiter.join(self.colname) + "\n"
        outfile = StringIO()
        np.savetxt(outfile, self.data, fmt, delimiter)
        output += outfile.getvalue()
        outfile.close()
        if FILENAME == None: return output
        else: 
            fh = open(FILENAME, "w")
            fh.write(output)
            fh.close()

