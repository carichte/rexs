from core import *
import numpy as np
import evaluationtools as et
import deltaf
from scipy import interpolate

def transform(element, energy, f2=None, f1=None, edge=None):
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
        raise ValueError("No fine structure given. Specify either `f1' or `f2'"
                         " keyword argument")
    
    
    dE = np.diff(energy)
    if np.allclose(dE, dE[0]):
        energy_org = energy
    else:
        energy_org = energy.copy()
        ffunc = interpolate.interp1d(energy, f)
        energy = np.arange(energy[0], energy[-1], dE.min())
        f = ffunc(energy)
    if f1:
        fsmooth, fTsmooth = deltaf.subtract_smooth(energy, element, f1=f,
                                                   edge=edge)
        ET, fT = imag(energy, f - fsmooth)
    elif f2:
        fsmooth, fTsmooth = deltaf.subtract_smooth(energy, element, f2=f,
                                                   edge=edge)
        ET, fT = real(energy, f - fsmooth)
    fT += fTsmooth(ET)
    
    fTfunc = interpolate.interp1d(ET, fT, bounds_error=False)
    
    fT = fTfunc(energy_org)
    return fT




