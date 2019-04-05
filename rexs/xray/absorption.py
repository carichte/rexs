#!/usr/bin/env python
#
#
# this is part of the rexs package written by Carsten Richter
# (carsten.richter@desy.de)
#

import os

import numpy as np
from scipy import interpolate, optimize

from . import interactions as xi
from . import deltaf
from .. import tools
from ..tools import wrap4leastsq


class Absorption(object):
    """
        A class to handle fluorescence data and to simulate fluorescence intensity
        or to extract the absorption coefficient from measured fluorescence data.
        
        Further, the absorption in bragg reflection regime can be calculated.
        
        The equations rely on the paper: C H Booth and F Bridges 2005 Phys. Scr. 2005 202
        
        
    """
    const = 10135467.657934014 # 2*eV/c/hbar
    def __init__(self, composition, resatom, energy, emission_energy, 
                       density, dE=0., Eedge=None, table="Sasaki", fwhm_eV=1):
        """
            Initialize instance to simulate fluorescence.
            The sample composition must be known and a resonant atom has to be assigned.
            
            Inputs: (* indicates mandatory inputs)
                *composition : string
                    sum formula of the sample material
                *resatom : string, int
                    symbol or number of the resonantly excited atom
                *energy : numpy.array
                    energy values - independent values for all methods
                    in eV
                *emission_energy : float
                    fluorescence line energy in eV
                *density : float
                    density of the sample compound in g/cm^3
                dE : float
                    edge shift compared to isolated atom in eV
                Eedge : float
                    Edge energy (eV). Will be fetched from database by default.
                table : database for x-ray optical constants
                
                
        """
        if isinstance(resatom, int):
            resatom = deltaf.elements.Symbols[resatom]
        assert(resatom in composition), "Resonant atom not found in composition."
        
        assert (np.diff(energy) > 0).all(), "Need monotonically increasing energy values"
        
        self.composition = composition
        self.density = density
        self.energy = energy
        self.resatom = resatom
        self.table=table
        self.dE = dE
        self.fwhm_eV = fwhm_eV

        self._weights = dict(Fluorescence = np.zeros(len(energy)),
                             Transmission = np.zeros(len(energy)),
                             Reflection   = np.zeros(len(energy)),
                             Absorption   = np.zeros(len(energy)))
        self.solve_data = "Fluorescence" # the best data for which abs coeff will be solved for
        self.IFluo = None
        self.IBragg = None
        self.DAFScalc = None
        self.Transmission = None
        self.callback = None
        
        
        self.beta_nonres = xi.get_optical_constants(density, composition, 
                                energy - dE, table=table, feff={resatom:0})[1]
        self.beta_tot = xi.get_optical_constants(density, composition, 
                                energy - dE, table=table)[1]

        self.mu_fluo = xi.get_optical_constants(density, composition, 
                                                emission_energy, 
                                                table=table)[1]
        
        if fwhm_eV:
            self.beta_nonres = tools.lorentzian_filter1d(self.beta_nonres, fwhm_eV)
            self.beta_tot = tools.lorentzian_filter1d(self.beta_tot, fwhm_eV)
            self.mu_fluo = tools.lorentzian_filter1d(self.mu_fluo, fwhm_eV)
        
        self.mu_nonres = self.beta_nonres * self.const * energy # nonresonant part of absorption
        self.mu_tot = self.beta_tot * self.const * energy # total absorption
        self.mu_fluo *= self.const * emission_energy
        self.mu_res_tab = self.mu_tot - self.mu_nonres # resonant part of absorption
        
        if Eedge is None:
            edges = deltaf.get_edge(resatom)
            edges_inv = dict([(v,k) for (k,v) in edges.items()])
            Eedges = edges.values()
            Eedges = filter(lambda E: (E > energy[0]) * (E < energy[-1]), Eedges)
            if len(Eedges)>1:
                print("Found multiple edges in energy range. "
                      "Using lowest (first) one.")
            Eedge = min(Eedges)
            
            self.Edge = Edge = edges_inv[Eedge]
            print("Using energy of %s edge at %.1f eV"%(Edge, Eedge))
        
        self.Eedge = Eedge
        
        ind = (energy - dE) < (Eedge - 5)
        
        parabola = tools.PolynomialFit(np.log(energy), 
                                    np.log(self.mu_res_tab), 
                                    indf=ind)
        
        self.mu_res_tab -= np.exp(parabola) # only lowest shell absorption
        self.mu_nonres  += np.exp(parabola)
        
        self.mumax = self.mu_res_tab.max()
        self.mu = self.mu_res_tab.copy()#tools.lorentzian_filter1d(self.mu_res_tab, fwhm_eV)
        self.mu_tot = self.mu_nonres + self.mu
        self.muguess = self.mu.copy()
        # scan parameters:
        self.p = {"omega":0., "d":np.inf, "theta":45., "om_range":0,
                  "m":0, "n":1., "mf":0, "nf":1., "theta_fluo":45.,
                  "mt":0, "nt":1, "tdead":0}
        self.pname = {
            "omega":"miscut",
            "d":"sample thickness",
            "theta":"bragg angle",
            "theta_fluo":"angle of fluorescence detector",
            "om_range":"range of miscuts on wavy surface",
            "m":"slope of linear background in Bragg intensity",
            "n":"offset of linear background in Bragg intensity",
            "mf":"slope of linear background in fluorescence intensity",
            "nf":"offset of linear background in fluorescence intensity",
            "mt":"slope of linear background in transmission",
            "nt":"offset of linear background in transmission",
            "tdead":"Ratio of dead time to measurement time"
            }
        
        
        self.xval = np.linspace(-4, 4, 51)
        self.stdnorm = tools.standard_normal(self.xval)
    
    
    
    def getomega(self, omega=None, om_range=None):
        """
            Generates a gaussian distribution of omega values according to
            the om_range input argument.
            Returns a f(omega), omega tuple.
        """
        if omega is None:
            omega = self.p["omega"]
        if om_range is None:
            om_range = self.p["om_range"]
        if np.ndim(omega)==0:
            omega = np.array(omega, ndmin=2)
            if om_range!=0:
                omega = self.xval[:,np.newaxis] * om_range + omega
                parts = self.stdnorm[:,np.newaxis] * np.diff(self.xval)[0]
            else:
                parts = 1
        elif np.ndim(omega)==1:
            parts = np.ones(len(omega), dtype=float)[:,np.newaxis]
        elif np.ndim(omega)==2:
            parts = omega[:,0][:,np.newaxis]
            omega = omega[:,1][:,np.newaxis]
        else: 
            raise ValueError("omega has to many dimensions")
        parts/=np.sum(parts) # normalization
        
        return parts, omega
    
    
    
    def set_Reflection(self, Intensity, Icalc, weights=1):
        Intensity = np.array(Intensity, ndmin=1)
        Intensity /= Intensity.max()
        Icalc = np.array(Icalc, ndmin=1)
        Icalc /= Icalc.max()
        self.IBragg = Intensity
        self.DAFScalc = Icalc
        self._weights["Reflection"] = np.ones(len(self.energy)) * weights
    
    def set_Fluorescence(self, Intensity, weights=1, solve_it=True):
        Intensity = np.array(Intensity, ndmin=1)
        ind = (self.Eedge + 100)
        postedge = np.mean((self.Eedge, self.energy.max()))
        poly = tools.PolynomialFit(self.energy, Intensity / self.mu_res_tab, 
                               [postedge, self.energy.max()], 10)
        Intensity /= poly
        self.p["nf"] = np.median(Intensity)/np.median(self.calc_Fluorescence())
        self.IFluo = Intensity
        self._weights["Fluorescence"] = np.ones(len(self.energy)) * weights
        if solve_it:
            self.solve_data = "Fluorescence"

    def set_Transmission(self, Intensity, weights=1, solve_it=False):
        Intensity = np.array(Intensity, ndmin=1)
        Intensity /= Intensity.max()
        self.Transmission = Intensity
        self._weights["Transmission"] = np.ones(len(self.energy)) * weights
        if solve_it:
            self.solve_data = "Transmission"
    
    def set_mu_weights(self, weights=1):
        self._weights["Absorption"] = np.ones(len(self.energy)) * weights
    
    
    
    def calc_Fluorescence(self, mu_res=None, **kwargs):
        """
            Calculates relative Fluorescence from the given resonant part of
            the absorption coefficient (mu_res) and the geometry parameters in
            the self.p dictionary. 
            If not given, self.mu is taken per default. 
            
            keyword arguments:
                - all elements of self.p dictionary
        """
        self.p.update(dict([i for i in kwargs.items() if i[0] in self.p]))
        
        parts, omega = self.getomega()
        
        bgfluo = self.p["mf"] * (self.energy - self.energy[0]) \
                    + self.p["nf"] # lin background
        
        if mu_res is None:
            mu_res = self.mu
        
        mu_res = abs(mu_res)
        
        s_in = np.sin(np.radians(self.p["theta"] + omega))
        # sine of incoming beam angle
        if self.p.has_key("theta_fluo"):
            s_ou = np.sin(np.radians(self.p["theta_fluo"]))
        else:
            s_ou = np.sin(np.radians(self.p["theta"] - omega))
        
        Int = parts * mu_res \
              / (mu_res + self.mu_nonres + s_in/s_ou*self.mu_fluo) \
              * (1 - np.exp(-((self.mu_nonres + mu_res)/s_in \
                             + self.mu_fluo/s_ou) * self.p["d"]))
        Int = Int.sum(0)
        
        if self.p["tdead"] != 0:
            Int = 1./(1./Int + abs(self.p["tdead"]))
        return Int * bgfluo
    
    
    
    def calc_Reflection(self, mu_tot=None, pol="sigma", **kwargs):
        """
            Calculates absorption in Bragg geometry.
            
            All necessary inputs are taken from the Absorption instance.
            Optionally they can be overwritten by passing to the function:
                
                mu_tot : np.array
                       total absorption coefficient of the material
                
                pol :  polarization of the x-rays
                
                The kwargs are all content of the self.p geometry parameters.
        """
        self.p.update(dict([i for i in kwargs.items() if i[0] in self.p]))
        
        if pol=="sigma":
            LP = 1.
        else:
            pass #not yet
        parts, omega = self.getomega()
        
        # lin background:
        bg = self.p["m"] * (self.energy - self.energy[0]) + self.p["n"] 
        
        Q = LP
        if mu_tot is None:
            #mu_tot = self.solve_mu_tot()
            mu_tot = self.mu_tot
        t_om = np.tan(np.radians(omega))
        t_th = np.tan(np.radians(self.p["theta"]))
        s_in = np.sin(np.radians(omega + self.p["theta"]))
        
        
        Int = parts * Q / (2*mu_tot) * (1 - t_om/t_th) \
            * (1 - np.exp(-2*mu_tot*self.p["d"]/s_in/(1 - t_om/t_th)))
        Int = Int.sum(0)
        Int *= bg / Int[0]
        self.Abs = Int
        return Int * self.DAFScalc
    
    
    
    def calc_Transmission(self, mu_tot=None, **kwargs):
        """
            Calculates absorption in Transmission geometry.
            
            All necessary inputs are taken from the self.p geometry parameters
            and the absorption coefficient self.mu.
            Optionally they can be overwritten by passing to the function:
                
                mu_tot : np.array
                       total absorption coefficient of the material
                
                The kwargs are all content of the self.p geometry parameters.
        """
        self.p.update(dict([i for i in kwargs.items() if i[0] in self.p]))
        
        parts, omega = self.getomega()
        
        # lin background:
        bgtrans = self.p["mt"] * (self.energy - self.energy[0])\
                     + self.p["nt"] 
        
        if mu_tot is None:
            #if self.solve_data != "Transmission":
            #    mu_tot = self.solve_mu_tot()
            mu_tot = self.mu_tot
            
        s_in = np.sin(np.radians(omega + self.p["theta"]))
        
        d = self.p["d"] / s_in
        Int = np.exp(-d*mu_tot) * parts
        Int = Int.sum(0)
        Int *= bgtrans
        return Int
    
    
    
    def solve_mu_tot(self, muguess=None, **kwargs):
        """
            Solve equation for Fluorescence by point-wise variation of the 
            absorption coefficient for given geometry parameters in self.p to 
            meet the given Fluorescence or Transmission. 
            
            Returns total absorption coefficient.
        """
        self.p.update(dict([i for i in kwargs.items() if i[0] in self.p]))
        
        if muguess is None:
            muguess = self.muguess
            
        
        if self.solve_data == "Fluorescence":
            resfunc = lambda x: (self.calc_Fluorescence(x) - self.IFluo)**2
            self.mu = abs(optimize.fsolve(resfunc, muguess))
            self.mu_tot = self.mu + self.mu_nonres
        elif self.solve_data == "Transmission":
            if self.p["om_range"] > 0:
                raise ValueError("om_range not supported for Transmission")
            parts, omega = self.getomega()
            omega = omega.item()
            bgtrans = self.p["mt"] * (self.energy - self.energy[0]) \
                                   + self.p["nt"]
            s_in = np.sin(np.radians(omega + self.p["theta"]))
            d = self.p["d"] / s_in
            I = self.Transmission / bgtrans
            self.mu_tot = -np.log(I)/d
            self.mu = self.mu_tot - self.mu_nonres
        return self.mu_tot
    
    
    
    def getf2(self, mu_tot=None):
        """
            Returns imaginary part of scattering amplitude ``f'' for resonant
            atom and given composition.
        """
        if mu_tot is None:
            mu_tot = self.mu_tot
        
        beta_tot = mu_tot / self.const / self.energy
        elements, amount = xi.get_components(self.composition)
        atweights = [xi.get_element(element)[1]/1000. for element in elements]
        beta_tot *= (np.array(amount)*np.array(atweights)).sum()
        beta_tot /= xi.electron_radius/(2*np.pi) \
                    * (xi.keV_A*1e-7/self.energy)**2 \
                    * self.density*1000. * xi.avogadro
        for i in range(len(elements)):
            if elements[i]==self.resatom:
                ires = i
                continue
            f1, f2 = np.array(xi.get_f1f2_from_db(elements[i], 
                                                  self.energy - self.dE, 
                                                  table=self.table))
            if self.fwhm_eV:
                f1 = tools.lorentzian_filter1d(f1, self.fwhm_eV)
                f2 = tools.lorentzian_filter1d(f2, self.fwhm_eV)
            beta_tot  -= f2 * amount[i]
        
        self.f2 = beta_tot / amount[ires]
        return self.f2
    
    
    
    def residuals(self, **p):
        self.p.update(dict([i for i in p.items() if i[0] in self.p]))
        
        for key in ["d"]:
            self.p[key] = abs(self.p[key])
        
        
        res = []
        if self.callback is not None:
            self.callback()
        self.solve_mu_tot()############################################################################?
        
        
        w = self._weights["Reflection"]
        if self.IBragg is not None and self.DAFScalc is not None and w.sum() > 0.:
            res.append((self.IBragg - self.calc_Reflection()) * w)
            if self._relative_fit:
                res[-1] /= self.IBragg
        
        w = self._weights["Absorption"]
        if w.sum() > 0.:
            res.append((self.mu_res_tab - self.mu) * w / self.mumax)
            if self._relative_fit:
                ind = self.energy > self.Eedge
                res[-1][ind] /= (self.mu_res_tab / self.mumax)[ind]
        
        w = self._weights["Transmission"]
        if self.Transmission is not None and w.sum()>0.:
            res.append((self.Transmission - self.calc_Transmission()) * w)
            if self._relative_fit:
                res[-1] /= self.Transmission
        
        w = self._weights["Fluorescence"]
        if self.IFluo is not None and w.sum()>0.:
            res.append((self.IFluo - self.calc_Fluorescence()) * w)
            if self._relative_fit:
                ind = self.energy > self.Eedge
                res[-1][ind] /= self.IFluo[ind]
        
        res = np.hstack(res)
        err = (res**2).sum()
        if err<self.errmin:
            self.muguess = self.mu
            self.errmin = err
        self.err = err
        print((self.err, "\t", " ".join(["%s=%g"%i for i in self.p.items()])))
        if self.fitalg == "simplex":
            return err
        else:
            return res
        
    def fitit(self, variables, fitalg="leastsq", callback=None, relative=True, **kwargs):
        """
            Fits the calculated DAFS (self.DAFScalc) to the measured Bragg 
            intensity (IBragg) by varying the geometry parameters in 
            self.p
        """
        self._relative_fit = relative
        self.callback = callback
        self.fitalg = fitalg
        self.variables = variables
        self.err = np.inf
        self.errmin = np.inf
        func, startval = wrap4leastsq.wrap_for_fit(self.residuals, self.p, variables)
        if self.fitalg == "simplex":
            output = optimize.fmin(func, startval, full_output=True, **kwargs)
                   #, maxfun=1000*len(startval), maxiter=1000*len(startval))
        else:
            output = optimize.leastsq(func, startval, full_output=True, **kwargs)
        param = output[0]
        return param
        

mu2fluo = Absorption # Backwards compatability
