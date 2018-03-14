#!/usr/bin/env python
#
#
# this is part of the rexs package written by Carsten Richter
# (carsten.richter@desy.de)
#

import os
import numpy as np
from scipy import ndimage, optimize, special

from ..io import FIOdata
from ..tools import diverse
from .interactions import get_components
from . import deltaf

class SDDspectrum(object):
    lambda_c = 2.4263106000088493e-2 # Comptonwellenlaenge in Angstrom
    def __init__(self, composition=None, erange=None):
        """
            Class to model histograms obtained by energy dispersive x-ray
            detectors (via multi channel analyzers ``mca``).
            Fluorescence lines can be modeled with the asymmetric detector
            response function of silicon drift diodes which goes beyond
            the simple gaussian peak shape.
            
            Lines originating of the same edge are linked via fixed ratios
            and their intensities varied accordingly to fit the measured 
            spectra.
            
            
            Optional inputs:
                composition : str
                    String containg all element symbols that are interesting
                    for data analysis. E.g. a sum formula.
                
                erange : 2-tuple of floats
                    2 values definingt lower and upper limit of the energy
                    range of interest (in eV).
        """
        self.fano = 49.686 # FANO limit des Rauschens bei 5895eV
        self.bg = 0.
        self.preedge = 15. # width of preedge range left of edge
        self.elastic = None
        
        # Standardwerte?
        self.p = {'A': 3e-4, 'B': -10., 'K0': 0., 'X': 0., 'c': 1.,
                  'b': 1e-2, 'w': 36.}
        self.edges={}
        self.lines={}
        self.strengths={"elastic":1., "compton":0.2}
        if erange is not None and len(erange)==2:
            self.emin, self.emax = min(erange), max(erange)
        else:
            self.emin, self.emax = 0, np.inf
        if composition is not None:
            elements = get_components(composition)[0]
            for element in elements:
                #print element
                transitions = deltaf.get_transition(element,  col='Direct')
                pick = filter(lambda x: np.isfinite(transitions[x]), transitions.keys())
                pick = filter(lambda x: transitions[x]>self.emin and transitions[x]<self.emax, pick)
                edges = filter(lambda x: "edge" in x, pick)
                lines = filter(lambda x: "edge" not in x, pick)
                for edge in edges:
                    thislines = np.array(filter(lambda s: s.startswith(edge.split()[0]), lines))
                    thislineseV, ind = np.unique([transitions[line] for line in thislines], return_index=True)
                    thislines = thislines[ind]
                    thislines = map(lambda s: element + "_" + s, thislines)
                    self.add_edge(transitions[edge], element + "_" + edge.split()[0], thislineseV, names=thislines)
            
        
        
    def add_edge(self, Eedge, name, lines, names=None, strength=1.):
        """
            - Eedge: energy of absorption edge in eV
            - lines: list of energies of the resulting emission lines in eV
            - names: list of  names   of the resulting emission lines in eV (optional)
        """
        if not hasattr(lines, "__iter__"):
            lines = (lines,)
        if names is not None:
            if not hasattr(names, "__iter__"):
                names = (names,)
            assert len(names)==len(lines), 'lines and names are not of same length'
        else:
            names = ["Line_%s%.2i"%(name, i+1) for i in range(len(lines))]
        ratios = list(np.ones(len(lines)))
        
        self.edges[name] = Eedge
        self.lines[name] = map(lambda *x: list(x), names, lines, ratios)
        #self.lines[name] = zip(names, lines, ratios)
        self.strengths[name] = strength
    
    def detector_response(self, channels, E0, amp):
        w = self.p["c"] * np.sqrt(self.p["w"]**2 + self.fano**2 / 5895. * E0)
        C0 = self.p["c"] * E0 + self.p["K0"]
        try:
            return amp * self.p["A"]/2. * np.exp(-channels*self.p["b"]/self.p["B"])\
                                        * (special.erf((channels-C0-self.p["X"])/self.p["B"])+1)\
                                        + diverse.gaussian(channels, C0, abs(amp), w, 0)
        except:
            print C0, w, amp
            raise
    
    def __call__(self, x=None):
        if x is None: x = self.channels
        output = self.bg * np.ones(len(x))
        for edge in self.edges.keys():
            if self.edges[edge] < (self.elastic - self.preedge):
                for line in self.lines[edge]:
                    output += self.detector_response(x, line[1], self.strengths[edge]*line[2])
        if not self.elastic is None:
            output += self.detector_response(x, self.elastic, self.strengths["elastic"])
            if hasattr(self, "twotheta"):
                output += self.detector_response(x, self.Ecompton, self.strengths["elastic"]*self.strengths["compton"])
        return output
    
    def set_variables(self, shape=[], strengths=[], ratios=[]):
        """
            Sets the set of variables for fitting.
            Inputs:
                - shape: list, empty or contains strings of shape parameters for
                               calculation the detector response
                               (either of "A", "b", "B", "c", "K0", "w", "X")
                               (see self.p.keys())
                - strengths: list, empty or contains names of edges whose absorption
                                   strength shall be varied
                                   (see self.strengths.keys())
                - ratios: list, empty or contains names of edges whose emission
                                line intensity ratios shall be varied
                                (see self.lines.keys())
        """
        if shape=="all": 
            shape = self.p.keys()
            shape.pop(shape.index("X"))
        if strengths=="all": strengths = self.strengths.keys()
        if ratios=="all": ratios = self.lines.keys()
        self.ShapeVars = []
        self.StrengthVars = []
        self.RatioVars = []
        self.guess = []
        for key in shape:
            if key in self.p.keys():
                self.ShapeVars.append(key)
                self.guess.append(self.p[key])
        for key in strengths:
            if self.strengths.has_key(key):
                if (key in self.edges.keys() and self.edges[key] < (self.elastic - self.preedge))\
                or  key in ["compton", "elastic"]:
                    self.StrengthVars.append(key)
                    self.guess.append(self.strengths[key])
            else:
                print "Warning: Edge `%s` is not defined!" %key
        for key in ratios:
            if self.edges.has_key(key)\
            and self.edges[key] < (self.elastic - self.preedge):
                self.RatioVars.append(self.lines[key])
                # ratio of first line always 1:
                self.guess.extend([line[2] for line in self.lines[key][1:]])
            else:
                print "Warning: Edge `%s` is not defined!" %key
    
    def fitfunction(self, *t):
        x = t[0]
        params = t[1:]
        i = len(self.ShapeVars)
        j = len(self.StrengthVars)
        self.p.update(dict(zip(self.ShapeVars, params[0:i])))
        self.strengths.update(dict(zip(self.StrengthVars, params[i:(i+j)])))
        for edge in self.RatioVars:
            for line in edge[1:]: # ratio of first line always 1
                line[2] = params[j]
                j+=1
        return self.__call__(x)
    
    def parse_mca(self, data, twotheta, energy=None, col=1, chmin=0, chmax=np.inf,
                        average=0, verbose=False):
        """
            Loads and processes .fio files produced by the ``online``
            software at DESY-FS.
            Here the data is supposed to be 1-dimensional, meaning 
            that only one COLUMN of the file will be processed.
            
            Inputs:
                data : str or file handle
                    Name or file handle of the .fio file to load
                
                twotheta : str or float
                    Either value or motor name containing the angle 
                    between incident and scattered beam.
                
                col : int
                    Number of column that shall be processed if data 
                    is 2-dimensional
                
                chmin, chmax : int
                    Channel limits whithin that the data will be 
                    cropped.
                
                average : int
                    Width of box filter for smoothing data in channels
                
                verbose : bool
                    Talk a lot?
        """
        if isinstance(data, np.ndarray):
            if data.ndim == 1:
                self.mcadata = data.copy()
            elif data.ndim == 2:
                self.mcadata = data[:,col]
        else:
            self.mca= FIOdata(data, verbose=verbose)
            self.mcadata = self.mca[:,col]
        self.channels = np.arange(len(self.mcadata)).astype(float)
        self.ind = (self.channels>chmin) * (self.channels<chmax)
        ind2 = self.mcadata < self.bg
        self.mcadata[ind2] = self.bg
        try:
            self.Iint = self.mcadata.sum() / self.mca.parameters["SAMPLE_TIME"]
            #print "Integral Intensity = %2f cps" %self.Iint
        except:
            self.Iint = self.mcadata.sum()
            #print "Integral Intensity = %2f counts" %self.Iint
        self.energy = (self.channels - self.p["K0"]) / self.p["c"]
        if average:
            self.mcadata = ndimage.uniform_filter1d(self.mcadata.astype(float), average)
        if energy is None and hasattr(self, "mca"):
            try:
                self.elastic = self.mca.parameters["ENERGY"]
            except KeyError:
                print("Warning: energy of incident photons could not be retrieved.")
        else:
            self.elastic = energy
        try: self.twotheta = float(twotheta)
        except: self.twotheta = self.mca.parameters[twotheta]
        Lcompton = self.lambda_c * (1 - np.cos(np.radians(self.twotheta))) \
                   + 12398./self.elastic
        self.Ecompton = 12398. / Lcompton # Energy of compton peak
    
    def refresh(self):
        self.energy = (self.channels - self.p["K0"]) / self.p["c"]
    
    def feed_energy(self, energy):
        energy = np.array(energy).ravel()
        if not hasattr(self, "channels"):
            self.channels = np.arange(len(energy))
        assert energy.shape == self.channels.shape, \
           "Given energy spectrum has to have the same length as given mca spectrum"
        self.energy = energy.copy()
        # channels = energy * c + K0
        c = (self.channels[-1] - self.channels[0]) / (energy[-1] - energy[0])
        K0 = self.channels[0] - c * energy[0]
        self.p.update(dict(c=c, K0=K0))
        
    
    def do_fit(self, verbose=False):
        assert hasattr(self, "guess"), "No set of variables defined. Run self.set_variables() first."
        assert hasattr(self, "mcadata"), "No measured dataset loaded. Run self.parse_mca() first."
        ResidualFunction = lambda t: (self.mcadata[self.ind] \
                                      - self.fitfunction(self.channels[self.ind], *t)) \
                                     / np.sqrt(self.mcadata[self.ind]+1)
        output = optimize.leastsq(ResidualFunction, self.guess, full_output=True, ftol=2**-20, xtol=2**-20)
        if np.isscalar(output[0]): fp = [output[0]]
        else: fp = output[0]
        for i in range(len(self.ShapeVars), len(self.StrengthVars) + len(self.RatioVars)):
            fp[i] = abs(fp[i])
        err = (ResidualFunction(fp)**2).sum()
        if output[1] is None: stddev = [np.inf for i in range(len(self.guess))]
        else: stddev = [np.sqrt(var) for var in output[1].diagonal()] # Kovarianzmatrix
        if verbose: print "Error at minimum: %f" %err
        return fp, stddev, err

