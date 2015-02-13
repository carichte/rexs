#!/usr/bin/env python
#
#
# this is part of the 'evaluationtools' set written by Carsten Richter
# (carsten.richter@desy.de)
#

import os
import evaluationtools as et
import pyxrr.xray_interactions as xi
#import mskk
import numpy as np
from wrap4leastsq import wrap_for_fit
from scipy import interpolate, optimize
import deltaf


class Absorption(object):
    """
        A class to handle fluorescence data and to simulate fluorescence intensity
        or to extract the absorption coefficient from measured fluorescence data.
        
        Further, the absorption in bragg reflection regime can be calculated.
        
        The equations rely on the paper: C H Booth and F Bridges 2005 Phys. Scr. 2005 202
        
        
    """
    const = 10135467.657934014 # 2*eV/c/hbar
    def __init__(self, composition, resatom, energy, emission_energy, 
                       density, dE=0., Eedge=None, table="Sasaki"):
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
        
        self.mu_nonres = self.beta_nonres * self.const * energy # nonresonant part of absorption
        self.mu_tot = self.beta_tot * self.const * energy # total absorption
        self.mu_fluo = xi.get_optical_constants(density, composition, 
                                                emission_energy, 
                                                table=table)[1]
        self.mu_fluo *= self.const * emission_energy
        self.mu_res_tab = self.mu_tot - self.mu_nonres # resonant part of absorption
        
        if Eedge == None:
            edges = deltaf.get_edge(resatom)
            edges_inv = dict([(v,k) for (k,v) in edges.iteritems()])
            Eedges = edges.values()
            Eedges = filter(lambda E: (E > energy[0]) * (E < energy[-1]), Eedges)
            if len(Eedges)>1:
                print("Found multiple edges in energy range. "
                      "Using lowest (first) one.")
            Eedge = min(Eedges)
            
            self.Edge = Edge = edges_inv[Eedge]
            print("Using energy of %s edge at %.1f eV"%(Edge, Eedge))
        
        self.Eedge = Eedge
        
        ind = energy < (Eedge - 5)
        
        parabola = et.PolynomialFit(np.log(energy), 
                                    np.log(self.mu_res_tab), 
                                    indf=ind)
        
        self.mu_res_tab -= np.exp(parabola) # only lowest shell absorption
        self.mu_nonres  += np.exp(parabola)
        
        self.mumax = self.mu_res_tab.max()
        self.mu = et.lorentzian_filter1d(self.mu_res_tab, 1)
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
        self.stdnorm = et.standard_normal(self.xval)
    
    
    
    def getomega(self, omega=None, om_range=None):
        """
            Generates a gaussian distribution of omega values according to
            the om_range input argument.
            Returns a f(omega), omega tuple.
        """
        if omega==None:
            omega = self.p["omega"]
        if om_range==None:
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
        poly = et.PolynomialFit(self.energy, Intensity / self.mu_res_tab, 
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
        self.p.update(dict([i for i in kwargs.iteritems() if i[0] in self.p]))
        
        parts, omega = self.getomega()
        
        bgfluo = self.p["mf"] * (self.energy - self.energy[0]) \
                    + self.p["nf"] # lin background
        
        if mu_res == None:
            mu_res = self.mu
        if mu_res == None:
            #if self.solve_data != "Fluorescence":
            #self.solve_mu_tot()
            self.mu_tot
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
        self.p.update(dict([i for i in kwargs.iteritems() if i[0] in self.p]))
        
        if pol=="sigma":
            LP = 1.
        else:
            pass #not yet
        parts, omega = self.getomega()
        
        # lin background:
        bg = self.p["m"] * (self.energy - self.energy[0]) + self.p["n"] 
        
        Q = LP
        if mu_tot==None:
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
        self.p.update(dict([i for i in kwargs.iteritems() if i[0] in self.p]))
        
        parts, omega = self.getomega()
        
        # lin background:
        bgtrans = self.p["mt"] * (self.energy - self.energy[0])\
                     + self.p["nt"] 
        
        if mu_tot==None:
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
        self.p.update(dict([i for i in kwargs.iteritems() if i[0] in self.p]))
        
        if muguess == None:
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
        if mu_tot==None:
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
            beta_tot  -= f2 * amount[i]
        
        self.f2 = beta_tot / amount[ires]
        return self.f2
    
    
    
    def residuals(self, **p):
        self.p.update(dict([i for i in p.iteritems() if i[0] in self.p]))
        
        for key in ["d"]:
            self.p[key] = abs(self.p[key])
        
        
        res = []
        if self.callback!=None:
            self.callback()
        self.solve_mu_tot()############################################################################
        
        w = self._weights["Reflection"]
        if self.IBragg != None and self.DAFScalc!=None and w.sum() > 0.:
            res.append((self.IBragg - self.calc_Reflection()) * w)
        
        w = self._weights["Absorption"]
        if w.sum() > 0.:
            res.append((self.mu_res_tab - self.mu) * w / self.mumax)
        
        w = self._weights["Transmission"]
        if self.Transmission!=None and w.sum()>0.:
            res.append((self.Transmission - self.calc_Transmission()) * w)
        
        w = self._weights["Fluorescence"]
        if self.IFluo!=None and w.sum()>0.:
            res.append((self.IFluo - self.calc_Fluorescence()) * w)
        
        res = np.hstack(res)
        err = (res**2).sum()
        if err<self.errmin:
            self.muguess = self.mu
            self.errmin = err
        self.err = err
        print self.err, "\t", " ".join(["%s=%g"%i for i in self.p.iteritems()])
        if self.fitalg == "simplex":
            return err
        else:
            return res
        
    def fitit(self, variables, fitalg="leastsq", callback=None, **kwargs):
        """
            Fits the calculated DAFS (self.DAFScalc) to the measured Bragg 
            intensity (IBragg) by varying the geometry parameters in 
            self.p
        """
        self.callback = callback
        self.fitalg = fitalg
        self.variables = variables
        self.err = np.inf
        self.errmin = np.inf
        func, startval = wrap_for_fit(self.residuals, self.p, variables)
        if self.fitalg == "simplex":
            
            output = optimize.fmin(func, startval, full_output=True, **kwargs)
                   #, maxfun=1000*len(startval), maxiter=1000*len(startval))
        else:
            output = optimize.leastsq(func, startval, full_output=True, **kwargs)
        param = output[0]
        return param
        

mu2fluo = Absorption # Backwards compatability

#
#def fluo_to_mu(composition, density, energy, fluorescence, order=1, dE=None,
#               full_output=True, table="Sasaki"):
#    """
#        DEPRECATED--
#        Fits a given fluorescence curve to the tabulated absorption 
#        coefficient of any material for a given energy in the x-ray regime
#        given in eV.
#        Since the relation abscoeff(fluo) is rather non-trivial a polynomial 
#        dependence up to a certain order<=2 is assumed. Furthermore, an energy
#        dependence up to the same order is assumed du take account for the 
#        energy dependent device function (absorption, source performance etc.).
#        
#        Inputs:
#            composition : string
#                Sum formula of the material.
#            
#            density : float
#                Density of the material.
#            
#            energy : numpy.array
#                Array of energy values.
#            
#            fluorescence : numpy.array
#                Array of fluorescence intensities.
#            
#        Optional inputs:
#            - order : int
#                Order of the polynomial describing device function and 
#                fluorescence vs. absorption dependence
#            
#            - dE : float
#                Shift of edge position. Difference between measured and
#                theoretical edge position.
#            
#            - full_output : bool
#                If true, a dictionary containing additional information will
#                be returned as second return value.
#            
#            - table : string
#                Database table to be used for smooth dispersion corrections
#                of free atom. See pyxrr.functions.get_optical_constants.__doc__
#                for details.
#        
#        Returns:
#            - mu : the fitted curve
#            
#            if full_output == True further output is generated:
#                - dE : energy shift compared to the tabulated curves
#                    (should be some few eV)
#                    
#                - mud_guess : the curve generated by the inital parameters.
#                
#                - mu_tab : the smooth absorption coefficient taken from the.
#                    database.
#                
#                - param : the parameters that minimize the error function.
#                
#                - error : a list containing the value of the error function
#                    for each iteration.
#                
#                - weights : weighting automatically generated by use of the
#                    edge position.
#    """
#    fitfunc = optimize.fmin
#    maxiter = 10000
#    const = 10135467.657934014 # 2*e/c/hbar
#    energy_ext = np.arange(energy[0]-100, energy[-1]+100, np.diff(energy).min())
#    delta, beta = xi.get_optical_constants(density, composition, energy_ext, 
#                                           table=table)
#    mu = et.lorentzian_filter1d(const * energy_ext * beta, 2)
#    #fluorescence = et.lorentzian_filter1d(fluorescence, 1)
#    
#    fluorescence -= fluorescence.min()
#    
#    dmu = np.diff(const * energy_ext * beta)
#    ind = abs(dmu) > 10*abs(np.median(dmu))
#    Eedge = energy_ext[np.where(ind==True)]
#    ind = np.zeros(len(energy), dtype=bool)
#    first = 1
#    for E in Eedge:
#        ind += (energy>(E - 20.)) * (energy<(E + 100.))
#        if first:
#            if dE==None:
#                imax = np.diff(fluorescence[ind]).argmax()
#                dE =  energy[ind][[imax, imax+1]].mean() - E
#            ind = (energy>(E + dE - 20.)) * (energy<(E + dE + 100.))
#            first = 0
#            print("Energy shift: %.1f eV"%dE)
#        print("Edge at %.2f eV"%E)
#    weights = np.ones(len(energy))
#    weights[~ind] *= 15.
#    mufunc = interpolate.interp1d(energy_ext + dE, mu, bounds_error=False, 
#                                  fill_value=np.nan)
#    
#    def mud(a0, a1, b0, ab0, b1, ab1):
#        #return (1 + b1*energy + c0*energy**2) * \
#        #       (a0 + b2*fluorescence + c1*fluorescence**2 + d1*fluorescence**3)
#        return a0 \
#             + a1 * fluorescence \
#             + b0 * (energy - energy[0]) \
#             + ab0 * (energy - energy[0]) * fluorescence \
#             + b1 * (energy - energy[0])**2 \
#             + ab1 * (energy - energy[0])**2 * fluorescence
#    
#    errs = []
#    def residuals(*v):
#        if hasattr(v[0], "__iter__"): v = np.append(v[0], v[1:]).ravel()
#        res = (mufunc(energy) -  mud(*v)) * weights
#        errs.append((res**2).sum())#, v #show error on every step
#        #return res
#        return (res**2).sum()
#    
#    #fp = optimize.fmin(residuals, [0,0.,1.,0.])[0]
#    
#    guess = [mu.min(),
#             (mu.mean() - mu.min())/fluorescence.mean(),
#             np.median(dmu)/np.median(np.diff(energy)),
#             0.]
#    
#    if order>=1:
#        fp = fitfunc(residuals, 
#                     guess,
#                     args = [ 0., 0.], maxiter=maxiter, maxfun=maxiter)
#        fp = np.append(fp, [0,0])
#    if order>=2:
#        #fp = np.append(guess, [0,0])
#        fp = fitfunc(residuals, fp, maxiter=maxiter, maxfun=maxiter)
#    #print "Fit result:\n", fp
#    if full_output:
#        fullout = {"mu_tab": mufunc(energy), 
#                   "mud_guess": mud(*(guess + [0,0])), 
#                   "dE":dE,
#                   "param":fp, 
#                   "error":np.array(errs), 
#                   "weights":weights}
#        return mud(*fp), fullout
#    else:
#        return mud(*fp)
#

def rebin_k_space(energy, xafs, edge, dist=10, conststep=True):
    """
        A smoothing filter for EXAFS data. It works as a moving 
        box average with variable window size corresponding to a
        constant with in electron k-space.
        
        Inputs:
            energy : numpy.array
                Array of energy values.
            
            xafs : numpy.array
                Array of intensity values.
            
            edge : float
                Absorption edge energy. Here the electron wave
                vector k equals zero.
            
            dist : float
                distance beyond the edge, where smoothing shall
                start.
            
            conststep : bool
                True if returned result shall be equally spaced
                in energy instead of k
            
        Outputs:
            A 2-tuple containing the resulting, smoothed
            (energy, xafs) arrays.
        
    """
    newene = energy.copy()
    a = dist
    ind = energy > (edge + dist)
    newene[ind] = np.sqrt(2.*a*(energy[ind] - edge) - a**2) + edge
    
    data = np.vstack((newene, xafs)).T
    esmooth, xafssmooth = xi.rebin_data(data[ind], np.diff(newene[~ind])[0]).T
    esmooth = edge + ((esmooth - edge)**2 + a**2)/(2*a)
    esmooth = np.append(energy[~ind], esmooth)
    xafssmooth = np.append(xafs[~ind], xafssmooth)
    #ind = esmooth > (edge + dist)
    if conststep:
        xafsfunc = interpolate.interp1d(esmooth, xafssmooth, kind="cubic", bounds_error=False)
        #xafssmooth = xafsfunc(energy)
        xafs = xafsfunc(energy)
        print xafs.shape
        #xafs = np.append(xafs[~ind], xafssmooth)
        return energy, xafs
    else:
        return esmooth, xafssmooth


