import numpy as np
from scipy import interpolate
from ... import tools
from .. import deltaf
from .core import *

#def transformonce(element, energy, f2=None, f1=None, edge=None, corr=True):
#    """
#        Perform the Kramers Kronig Transform of the given real (f1) or 
#        imaginary (f2) part of the dispersion correction of the scattering 
#        amplitude `f'.
#    """
#    if f1 != None:
#        f = f1
#        f1 = True
#    elif f2 != None:
#        f = f2
#        f2 = True
#    else:
#        raise ValueError("No fine structure given. Specify either `f1' or `f2'"
#                         " keyword argument")
#    
#    
#    dE = np.diff(energy)
#    if np.allclose(dE, dE[0]):
#        energy_org = energy
#    else:
#        energy_org = energy.copy()
#        ffunc = interpolate.interp1d(energy, f)
#        energy = np.arange(energy[0], energy[-1], dE.min())
#        f = ffunc(energy)
#    if f1:
#        fsmooth, fTsmooth = deltaf.get_smooth(energy, element, f1=f,
#                                                   edge=edge)
#        ET, fT = imag(energy, f - fsmooth, corr=corr)
#    elif f2:
#        fsmooth, fTsmooth = deltaf.get_smooth(energy, element, f2=f,
#                                                   edge=edge)
#        ET, fT = real(energy, f - fsmooth, corr=corr)
#    fT += fTsmooth(ET)
#    
#    fTfunc = interpolate.interp1d(ET, fT, bounds_error=False)
#    
#    fT = fTfunc(energy_org)
#    return fT



class Transform(object):
    """
        Returns the smooth part either of the real (f1) or the imaginary
        (f2) part of the x-ray dispersion correction of the scattering
        amplitude f.
        The smooth part is negative for f1 and positive for f2 and 
        f1(E) -> 0 for E -> inf.
    """
    def __init__(self, energy, element, edges=(None,), fwhm=5., weights=None):
        self.energy = energy
        if not hasattr(edges, "__iter__"):
            edges = (edges,)
        Esmooth = []
        for i, edge in enumerate(edges):
            edgetab = deltaf.get_edge(element)
            if edge!=None and edge in edgetab:
                eedge = edgetab[edge]
            else:
                eedge = None
            smth = deltaf.get_energies(element, energy[0], energy[-1], 
                                       fwhm, eedge)
            Esmooth.append(smth[0])
            if i==0:
                Eedge = smth[0][smth[1]]
        
        self.Eedge = Eedge
        Esmooth = np.hstack(Esmooth)
        Esmooth.sort()
        self.Esmooth = Esmooth
        f1s = deltaf.getfquad(element, Esmooth, fwhm, f1f2 = "f1",
                              kernel="normal")
        self._f1func = interpolate.interp1d(Esmooth, f1s, kind="linear", 
                                            bounds_error=False, fill_value=0.)
        f2s = deltaf.getfquad(element, Esmooth, fwhm, f1f2 = "f2", 
                              kernel="normal")
        self._f2func = interpolate.interp1d(Esmooth, f2s, kind="linear", 
                                            bounds_error=False, fill_value=0.)
        if weights==None:
            weights = np.ones(len(energy))
            weights[5:-5]/= 10.
        self.weights=weights
        
    def fitfunc(self, E, func, E0, scale):
        return scale*func(E-E0)
        
    def get_smooth_parts(self, f1=None, f2=None):
        if f1==None and f2==None:
            print("No fine structure given. Specify either `f1' or `f2'"
                  " keyword argument")
            return None
        
        guess = dict(E0=0, scale=1)
        if f1 != None:
            E0m = self.energy[f1.argmin()]
            print("approximate measured edge: %g"%E0m)
            guess["E0"] = E0m - self.Eedge
            guess["func"] = self._f1func
            fit = tools.fitls(self.energy, f1, self.fitfunc, guess, 
                           ["scale", "E0"], fitalg="simplex", 
                           weights=self.weights)
            fit.popt["func"] = self._f2func
        elif f2 != None:
            print np.diff(f2).argmax()
            E0m = self.energy[np.diff(f2).argmax()]
            print("approximate measured edge: %g"%E0m)
            guess["E0"] = E0m - self.Eedge
            guess["func"] = self._f2func
            fit = tools.fitls(self.energy, f2, self.fitfunc, guess, 
                           ["scale", "E0"], fitalg="simplex", 
                           weights=self.weights)
            fit.popt["func"] = self._f1func
        print("Fit results:")
        print("scaling factor: %g"%fit.popt["scale"])
        print("edge shift: %g eV"%fit.popt["E0"])
        self.fit = fit
        return self.smoothf1, self.smoothf2
    
    def smoothf1(self, E):
        if not hasattr(self, "fit"):
            print("No smooth part found. Perform Transform.get_smooth_parts"
                  " first.")
            return None
        fp = self.fit.popt
        return self.fitfunc(E, self._f1func, fp["E0"], fp["scale"])
    
    def smoothf2(self, E):
        if not hasattr(self, "fit"):
            print("No smooth part found. Perform Transform.get_smooth_parts"
                  " first.")
            return None
        fp = self.fit.popt
        return self.fitfunc(E, self._f2func, fp["E0"], fp["scale"])
    
    def transform(self, f1=None, f2=None, corr=True):
        """
            Perform the Kramers Kronig Transform of the given real (f1) or 
            imaginary (f2) part of the dispersion correction of the scattering 
            amplitude `f'.
        """
        if f1 != None:
            f = f1
            f1 = True
        elif f2 != None:
            f = f2
            f2 = True
        else:
            raise ValueError("No fine structure given. Specify either `f1' "
                             "or `f2' keyword argument")
        
        if f1:
            fsmooth, fTsmooth = self.get_smooth_parts(f1=f)
        elif f2:
            fTsmooth, fsmooth = self.get_smooth_parts(f2=f)
        
        dE = np.diff(self.energy)
        if np.allclose(dE, dE[0]):
            energy = self.energy
        else:
            ffunc = interpolate.interp1d(self.energy, f)
            energy = np.arange(self.energy[0], self.energy[-1], dE.min())
            f = ffunc(energy)
        if f1:
            ET, fT = imag(energy, f - fsmooth(energy), corr=corr)
        elif f2:
            ET, fT = real(energy, f - fsmooth(energy), corr=corr)
        
        fT += fTsmooth(ET)
        
        fTfunc = interpolate.interp1d(ET, fT, bounds_error=False)
        fT = fTfunc(self.energy)
        return fT
    
