import numpy as np
import os
import elements
import deltaf
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.integrate import quad
from scipy.optimize import fmin
import evaluationtools as et

__doc__ =  """
        This is a set of routines for applications in resonant x-ray scattering.
        It includes interfaces to x-ray libraries and calculation of dispersion
        corrections.
    """

_std_fwhm_fac = 2 * np.sqrt(2*np.log(2))
_stwopi = np.sqrt(2*np.pi)

def cauchy(x, fwhm=1):
    w = fwhm/2.
    return 1./np.pi * w / (w**2 + x**2)

    
def normal(x, fwhm=1):
    sigma = fwhm / _std_fwhm_fac
    return np.exp(-(x/sigma)**2/2)/_stwopi/sigma

def get_transition(element, transition=None, col="Direct"):
    """
        Returns all x-ray absorption edge and line enegies for a
        given element. Optionally a certain transition can be chosen.
        
        Inputs:
            element : str
                exact symbol of the chemical element of interest.
        
        Optional inputs:
            transition : str
                name of selected transition
                possible values are:
                    ['K edge',
                     'KL1', 'KL2', 'KL3',
                     'KM1', 'KM2', 'KM3', 'KM4', 'KM5', 'KN1',
                     'L1 edge',
                     'L1M1', 'L1M2', 'L1M3', 'L1M4', 'L1N1',
                     'L2 edge', 'L2M1', 'L2M2', 'L2M3', 'L2M4', 'L2N1',
                     'L3 edge', 'L3M1', 'L3M2', 'L3M3', 'L3M4', 'L3M5', 'L3N1']
            
            col : str
                selected dataset.
                can be one of:
                    ["Theory", "TheoryUnc",
                     "Direct", "DirectUnc",
                     "Combined", "CombinedUnc",
                     "Vapor",  "VaporUnc"]
                where Unc denotes uncertainties of the values
        
        The data was obtained from (see for more information):
            Deslattes R D, Kessler E G, Indelicato P, de Billy L, Lindroth E and Anton J,
            Rev. Mod. Phys. 75, 35-99 (2003)
            10.1103/RevModPhys.75.35
            """
    DB_NAME = "transition_emission_energies.npz"
    DB_PATH = os.path.join(os.path.dirname(__file__), DB_NAME)
    try: Z = int(element)
    except: Z = elements.Z[element]
    db = np.load(DB_PATH)
    try:
        ind = db["Z"]==Z
        transitions = db["type"][ind]
        energies = db[col][ind]
    except Exception as errmsg:
        print errmsg
    finally:
        if hasattr(db, "close"):
            db.close()
    ind = np.where(transitions==transition)[0]
    if len(ind):
        return energies[ind]
    else:
        return dict(zip(transitions, energies))
    

def get_edge(element, edge=None):
    """
        Returns x-ray absorption edge of a chemical element of choice.
        A particular absorption edge can be chosen out of
            {Z, A, K, L1, L2, L3, M1, M5},
        where Z and A denote atomic number and atomic weight, respectively.
    """
    DB_NAME = "elements.npz"
    DB_PATH = os.path.join(os.path.dirname(__file__), DB_NAME)
    try: Z = int(element)
    except: Z = elements.Z[element]
    db = np.load(DB_PATH)
    cols = db.files
    try:
        ind = Z-1
        edges = dict([(col, db[col][ind]) for col in cols])
    except Exception as errmsg:
        print errmsg
    finally:
        db.close()
    if edge in cols:
        return edges[edge]
    else:
        return edges

def getf(element, energy, conv_fwhm=0):
    """
        This is a tiny wrapper for Matt Newville's deltaf Fortran program.
        It returns the correction terms f' and f'' of the atomic scattering
        factors for a given element and energy range.
        
        Inputs:
            - element : int or str - atomic number or abbreviation of element
            - energy  : array of energies in eV
            - conv_fwhm : float - width of lorentzian for broadening in eV.
        
        Outputs:
            f1, f2 Tuple
    """
    try: Z = int(element)
    except: Z = elements.Z[element]
    diffE = np.unique(np.diff(energy))
    dE = abs(diffE.min())
    fwhm_chan = abs(float(conv_fwhm))/dE/2.
    # make equally spaced data:
    if diffE.std()/diffE.mean() > 1e-5 and fwhm_chan>1e-2:
        newenergy = np.arange(energy.min() - 10*dE, energy.max() + 10.*dE, dE)
        do_interpolate = True
    else:
        newenergy = energy
        do_interpolate = False
        
    f1, f2 = deltaf.clcalc(Z,  newenergy)
    if fwhm_chan>1e-2:
        f1, f2 = deltaf.convl2(f1, f2, fwhm_chan)
    if do_interpolate:
        ffunc = interp1d(newenergy, (f1, f2))
        f1, f2 = ffunc(energy)
    return f1, f2
    

def getfquad(element, energy=None, fwhm_ev = 0, lw = 50, f1f2 = "f1",
             kernel="cauchy"):
    """
        This is a tiny wrapper for Matt Newville's deltaf Fortran program.
        It returns the correction terms f' or f'' of the atomic scattering
        factors for a given element and energy range.
        
        Here the energy doesn't need to be of equal steps since the integral
        is calculated numerically. This enables very fine stepping near absorption
        edges and coarse far from edges. 
        
        If energy is not given it is calculated on log-scale relative to the 
        edge. This way fine sampling is obtained with low number of points.
    """
    try: Z = int(element)
    except: Z = elements.Z[element]
    
    fwhm_ev = abs(fwhm_ev)
    
    if energy==None:
        energy, iedge = get_energies(Z, 1000, 10000, fwhm_ev=fwhm_ev)
        return_ene = True
    else:
        return_ene = False
    
    if f1f2 == "f2" or f1f2 == 1:
        ind = 1
    else:
        ind = 0
    
    
    if fwhm_ev <= np.finfo(float).eps:
        result =  deltaf.clcalc(Z, energy)[ind]
    else:
        
        if kernel=="cauchy":
            corrfac = 1./quad(cauchy, -lw, lw, args=(1,), limit=500)[0]
            integrand = lambda x, E: cauchy(x, fwhm_ev) * deltaf.clcalc(Z, E-x)[ind]
        elif kernel=="normal":
            corrfac = 1./quad(normal, -lw, lw, args=(1,), limit=500)[0]
            integrand = lambda x, E: normal(x, fwhm_ev) * deltaf.clcalc(Z, E-x)[ind]
        
        def ffunc(E):
            return quad(integrand, -lw*fwhm_ev, lw*fwhm_ev, args=(E,), limit=500)[0]
        
        if np.isscalar(energy) or len(energy)==1:
            result = corrfac * ffunc(energy)
        else:
            fvec = np.vectorize(ffunc)
            result = corrfac * fvec(energy)
    
    if return_ene:
        return energy, result
    else:
        return result



def get_energies(element, emin, emax, fwhm_ev=1e-4, eedge = None, num=100):
    """
        Returnes an array of energies in the vicinity of the given
        limits ``emin`` and ``emax`` so that the absorption edge of 
        ``element`` is sampled in a very fine way.
        The edge is found automatically.
        If there are many edges, the interesting edge can be given 
        by a rough estimate of the edge energy (``eedge``).
        
        A full width at half maximum value for energy resolution in eV can
        be given via the ``fwhm_ev`` keyword to conclude about the minimum 
        step width necessary. However it can not be zero, so the minimum
        is set internally to 1e-6 eV.
        
        The number of points on each side of the edge can be given by
        ``num``.
    """
    assert emax>emin, "emax must be larger than emin."
    fwhm_ev = max(abs(fwhm_ev), 1e-6)
    try: Z = int(element)
    except: Z = elements.Z[element]
    def f1func(E):
        return deltaf.clcalc(Z, E)[0]
    if eedge == None:
        eedge = .5 * (emin + emax)
    eedge = float(fmin(f1func, (eedge,)))
    print "Found edge at %g"%eedge
    expmin = np.floor(np.log10(fwhm_ev))
    dist = float(2*max(abs(emax - eedge), abs(emin - eedge))) * 1.2
    expmax = np.log10(dist/2.)
    expmin, expmax = min(expmin, expmax), max(expmin, expmax)
    ewing = np.logspace(expmin, expmax, num)
    dE = ewing[1] - ewing[0]
    num_lin = 2*int(ewing[0] / dE) + 1
    ecenter = np.linspace(-ewing[0], ewing[0], num_lin)[1:-1]
    energy = np.append(-ewing[::-1], ecenter)
    energy = np.append(energy, ewing)
    energy += eedge
    energy = energy[energy>0]
    try:
        iedge = float(np.where(energy==eedge)[0])
    except:
        raise ValueError("Edge energy not found in constructed array")
    return energy, iedge


def subtract_smooth(energy, element, f1=None, f2=None, fwhm=5., edge=None):
    """
        Subtracts the smooth part either of the real (f1) or the imaginary
        (f2) part of the x-ray dispersion correction of the scattering
        amplitude f.
        The smooth part is negative for f1 and positive for f2 and 
        f(E) -> 0 for E -> inf.
    """
    edges = get_edge(element)
    if edge in edges:
        eedge = edges[edge]
    else:
        eedge = None
    Esmooth, iedge = get_energies(element, energy[0], energy[-1], fwhm, eedge)
    eedge = Esmooth[iedge]
    
    f1s = getfquad(element, Esmooth, fwhm, f1f2 = "f1", kernel="normal")
    f1func = interp1d(Esmooth, f1s, kind="linear", bounds_error=False,
                      fill_value=0.)
    f2s = getfquad(element, Esmooth, fwhm, f1f2 = "f2", kernel="normal")
    f2func = interp1d(Esmooth, f2s, kind="linear", bounds_error=False,
                      fill_value=0.)
    
    
    def fitfunc(E, func, E0, scale):
        return scale*func(E-E0)
    
    
    
    width = min(abs(energy[0] - eedge), abs(energy[-1] - eedge))
    weights = np.ones(len(energy))
    weights[5:-5]/= 10.
    guess = dict(E0=0, scale=1)
    if f1 != None:
        E0m = energy[f1.argmin()]
        print("approximate measured edge: %g"%E0m)
        guess["E0"] = E0m - eedge
        guess["func"] = f1func
        fit = et.fitls(energy, f1, fitfunc, guess, ["scale", "E0"],
                       fitalg="simplex", weights=weights)
        fp = fit.popt
        fp["func"] = f2func
    elif f2 != None:
        E0m = energy[np.diff(f2).argmax()]
        print("approximate measured edge: %g"%E0m)
        guess["E0"] = E0m - eedge
        guess["func"] = f2func
        fit = et.fitls(energy, f2, fitfunc, guess, ["scale", "E0"],
                       fitalg="simplex", weights=weights)
        fp = fit.popt
        fp["func"] = f1func
    else:
        print("No fine structure given. Specify either `f1' or `f2' keyword"
              " argument")
    
        return None
    
    print("Fit results:")
    print("scaling factor: %g"%fit.popt["scale"])
    print("edge shift: %g eV"%fit.popt["E0"])
    return fit.ysim, fit.func


