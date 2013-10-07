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
from FIOdata import FIOdata
from diverse import *

class SDDspectrum(object):
    lambda_c = 2.4263106000088493e-2 # Comptonwellenlaenge in Angstrom
    def __init__(self):
        self.fano = 49.686 # FANO limit des Rauschens bei 5895eV
        self.bg = 0.
        self.elastic = None
        
        # Standardwerte?
        self.p = {'A': 3e-4, 'B': -10., 'K0': 0., 'X': 0., 'c': 1.,
                  'b': 1e-2, 'w': 36.}
        self.edges={}
        self.lines={}
        self.strengths={"elastic":1., "compton":0.2}
    def add_edge(self, edge, name, lines, names=None, strength=1.):
        """
            - edge: energy of absorption edge in eV
            - lines: list of energies of the resulting emission lines in eV
            - names: list of  names   of the resulting emission lines in eV (optional)
        """
        if not hasattr(lines, "__iter__"): lines = (lines,)
        if names is not None:
            if not hasattr(names, "__iter__"): names = (names,)
            assert len(names)==len(lines), 'lines and names are not of same length'
        else:
            names = ["Line_%s%.2i"%(name, i+1) for i in range(len(lines))]
        ratios = list(np.ones(len(lines)))
        
        self.edges[name] = edge
        self.lines[name] = map(lambda *x: list(x), names, lines, ratios)
        self.strengths[name] = strength
    
    def detector_response(self, channels, E0, amp):
        w = self.p["c"] * np.sqrt(self.p["w"]**2 + self.fano**2 / 5895. * E0)
        C0 = self.p["c"] * E0 + self.p["K0"]
        try:
            return amp * self.p["A"]/2. * np.exp(-channels*self.p["b"]/self.p["B"])\
                                        * (special.erf((channels-C0-self.p["X"])/self.p["B"])+1)\
                                        + gaussian(channels, C0, abs(amp), w, 0)
        except:
            print C0, w, amp
            raise
    
    def __call__(self, x=None):
        if x==None: x = self.channels
        output = self.bg * np.ones(len(x))
        for edge in self.edges.keys():
            if self.edges[edge] < self.elastic:
                for line in self.lines[edge]:
                    output += self.detector_response(x, line[1], self.strengths[edge]*line[2])
        if not self.elastic==None:
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
        if shape=="all": shape = self.p.keys()
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
            if key in self.strengths.keys():
                if (key in self.edges.keys() and self.edges[key] < self.elastic) or key in ["compton", "elastic"]:
                    self.StrengthVars.append(key)
                    self.guess.append(self.strengths[key])
            else: print "Warning: Edge `%s` is not defined!" %key
        for key in ratios:
            if (key in self.edges.keys() and self.edges[key] < self.elastic):
                self.RatioVars.append(self.lines[key])
                for line in self.lines[key][1:]: # ratio of first line alsways 1
                    self.guess.append(line[2])
    
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
    
    def parse_mca(self, fname, twotheta, average=10, col=1, chmin=0, chmax=np.inf, verbose=True):
        self.mca= FIOdata(fname, verbose=verbose)
        self.mcadata = self.mca[:,col]
        self.channels = (np.arange(len(self.mcadata)).astype(float)+1)
        self.ind = (self.channels>chmin) * (self.channels<chmax)
        ind2 = self.mcadata < self.bg
        self.mcadata[ind2]=self.bg
        try:
            self.Iint = self.mcadata.sum() / self.mca.parameters["SAMPLE_TIME"]
            #print "Integral Intensity = %2f cps" %self.Iint
        except:
            self.Iint = self.mcadata.sum()
            #print "Integral Intensity = %2f counts" %self.Iint
        self.energy = (self.channels - self.p["K0"]) / self.p["c"]
        self.mcadata = ndimage.uniform_filter1d(self.mcadata, average)
        self.elastic = self.mca.parameters["ENERGY"]
        try: self.twotheta = float(twotheta)
        except: self.twotheta = self.mca.parameters[twotheta]
        Lcompton = self.lambda_c * (1 - np.cos(np.radians(self.twotheta))) + 12398./self.elastic
        self.Ecompton = 12398. / Lcompton
        """
        for name in self.strengths.keys():
            if name=="elastic": 
                C0 = int(self.p["c"] * self.elastic + self.p["K0"])
                self.strengths["elastic"] = self.mcadata[C0]
                #self.strengths["compton"] = self.lines[
            elif name=="compton": continue
            else: 
                if self.edges[name]<self.elastic:
                    C0 = int(self.p["c"] * self.lines[name][0][1] + self.p["K0"])
                    self.strengths[name] = self.mcadata[C0]
                    for line in self.lines[name][1:]: # ratio of first line always 1
                        Ci = int(self.p["c"] * line[1] + self.p["K0"])
                        line[2] = self.mcadata[Ci] / self.mcadata[C0]
                else:
                    self.strengths[name] = self.mcadata.min()
        """
    
    def refresh(self):
        self.energy = (self.channels - self.p["K0"]) / self.p["c"]
    
    def do_fit(self, verbose=False):
        assert hasattr(self, "guess"), "No set of variables defined. Run self.set_variables() first."
        assert hasattr(self, "mcadata"), "No measured dataset loaded. Run self.parse_mca() first."
        ResidualFunction = lambda t: (self.mcadata[self.ind] \
                                      - self.fitfunction(self.channels[self.ind], *t)) \
                                     / np.sqrt(self.mcadata[self.ind])
        output = optimize.leastsq(ResidualFunction, self.guess, full_output=True, ftol=2**-20, xtol=2**-20)
        if np.isscalar(output[0]): fp = [output[0]]
        else: fp = output[0]
        for i in range(len(self.ShapeVars), len(self.StrengthVars) + len(self.RatioVars)):
            fp[i] = abs(fp[i])
        err = (ResidualFunction(fp)**2).sum()
        if output[1]==None: stddev = [np.inf for i in range(len(self.guess))]
        else: stddev = [np.sqrt(var) for var in output[1].diagonal()] # Kovarianzmatrix
        if verbose: print "Error at minimum: %f" %err
        return fp, stddev, err

