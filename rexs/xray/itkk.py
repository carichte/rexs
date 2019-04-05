import os
import numpy as np
import scipy.optimize as sopt
from scipy import interpolate, ndimage
from . import mskk
from . import deltaf

class IterativeKK(object):
    def __init__(self, cs, miller, energy, Idafs, f1start=None, f2start=None):
        """
            class for handling iterative Kramers-Kronig transformation which yields a
            Kramers-Kronig pair (f', f'') that describes the Diffraction Anomalous 
            Fine Strukture (DAFS) intensity (Idafs) measured in a certain energy range
            near an x-ray absorption edge.
            
            
            inputs:
                cs : the crystal structure given as an unit_cell instance
                
                miller : the bragg reflection of interest as 3-tuple
                
                energy : numpy array of photon energies in electron volt
                
                Idafs : measured and scaled DAFS intensity
        """
        self.cs = cs
        self.miller = miller
        self.energy = energy
        self.DAFS = Idafs
        for k in f1start:
            cs.feed_feff(k, energy, f1start[k], f2start[k])
        self.Ffunc = cs.DAFS(energy, miller, func_output=True)
        self.f, self.f1, self.f2 = {}, {}, {}
        self.f1tab, self.f2tab = {}, {}
        self.f1func, self.f2func = {}, {}
        self.Z = {}
        for symbol in cs.f:
            if symbol.name.startswith("f_"):
                self.f1[symbol] = cs.f[symbol].real.copy()
                self.f2[symbol] = cs.f[symbol].imag.copy()
                self.f[symbol.name] = cs.f[symbol].copy()
        #self.eneleft = 
        self.Isim = abs(self.Ffunc.dictcall(self.f))**2
        self.Isim0 = self.Isim.copy()
        self.ind={}
        self.ilim = {}
        self.anchors = {}
        self.diff = {}
        self.debug = True

    def AddFS(self, emin, emax, symbol, anchorpoints=None, w=10., f1func=None, f2func=None):
        element = symbol.name.split("_")[1]
        self.Z[symbol] = deltaf.elements.Z[element]
        #self.f1tab[symbol], self.f2tab[symbol] = xi.get_f1f2_from_db(element, self.energy - self.cs.dE[element], table="Sasaki")
        E, ind = deltaf.get_energies(element, self.energy[0], 
                                              self.energy[-1], fwhm_ev=w)
        eedge = E[ind]
        lw = min(max(abs(emax-eedge), abs(eedge-emin))*2./w, 200.)
        if f1func is None:
            #if os.path.isfile(".itkk_f1_%s.dmp")
            f1 = deltaf.getfquad(element, E, w, f1f2='f1', lw=lw) + self.Z[symbol]
            self.f1func[symbol] = interpolate.UnivariateSpline(E, f1, k=1, s=0)
        else:
            self.f1func[symbol] = f1func
        if f2func is None:
            f2 = deltaf.getfquad(element, E, w, f1f2='f2', lw=lw)
            self.f2func[symbol] = interpolate.UnivariateSpline(E, f2, k=1, s=0)
        else:
            self.f2func[symbol] = f2func
        self.diff[symbol] = complex(self.cs.F_0.subs(self.cs.subs).n().diff(symbol))
        ind = (self.energy > emin) * (self.energy < emax)
        self.ind[symbol] = ind
        
        ileft = abs(self.energy - emin).argmin()
        iright = abs(self.energy - emax).argmin()
        self.ilim[symbol] = (ileft, iright)
        
        if anchorpoints is None:
            pass
        else:
            f1f2func = interpolate.interp1d(self.energy, (self.f1[symbol], self.f2[symbol]), copy=False)
            ReAnch, ImAnch = f1f2func(anchorpoints)
            self.anchors[symbol] = (anchorpoints, ReAnch + 1j * ImAnch)
    
    def iterate(self, symbols=None, KK=None, overshoot=1):
        if KK is None:
            KK = {}
        E = self.energy
        dchan = 0.5/np.diff(E).mean() #0.5eV smear?
        if symbols is None:
            symbols = self.f1.iterkeys()
        allind = np.zeros(len(E), dtype=bool)
        for symbol in symbols:
            if symbol in self.anchors:
                ax,ay = self.anchors[symbol]
            else:
                ax, ay = None, None
            ind = self.ind[symbol]
            f1 = self.f1[symbol]
            f2 = self.f2[symbol]
            allind += ind
            facI = self.DAFS[ind] / self.Isim[ind]
            dF_F = np.sqrt(facI) - 1
            Fbefore = self.Ffunc.dictcall(self.f)[ind]
            F = (dF_F+1) * Fbefore
#            arg = ndimage.gaussian_filter1d(np.angle(F), dchan)
            arg = np.angle(F)
            F = abs(F)*np.exp(1j*arg)
            dF = F - Fbefore
            df = dF / self.diff[symbol] * overshoot
            self.f[symbol.name][ind] += df
            f1[ind] += df.real
            f2[ind] += df.imag
            self.F1 = self.Ffunc.dictcall(self.f)
            #self.F1 = F.copy()
            self.I1 = abs(self.F1)**2
            #err1 = ((self.DAFS[ind] - self.I1[ind])**2).sum()
            
            #if (abs(df.real) > abs(df.imag)).mean() > 0.5:
            if self.debug:
                import pylab as pl
                pl.plot(E, f1)
                pl.plot(E, f2)
            
            il, ir = self.ilim[symbol]
            """
            if abs(df.real).sum() > abs(df.imag).sum():
                KK[symbol] = "f1"
            else:
                KK[symbol] = "f2"
            """
            #if "Ti" in symbol.name or "Ba" in symbol.name:
            if KK.has_key(symbol) and KK[symbol]=="f1":
                KK[symbol] = "f1"
                if not ay is None:
                    ay = ay.imag
                """
                ancind = abs(E - np.array(ax)[:,np.newaxis]).argmin(-1)
                ay = f2[ancind]
                ax = E[ancind]
                """
                #thisfp = f1 - self.Z[symbol]
                thisf1 = f1 - self.f1func[symbol](E)
                newE, newf2 = mskk.imag(E, thisf1, ax, ay, corr=True)
                if self.debug:
                    import pylab as pl
                    pl.plot(E, f1, label="f1")
                    pl.plot(E, self.f1func[symbol](E), "--", label="f1tab")
                    pl.plot(E, thisf1, label="f1 - f1tab")
                    pl.plot(newE, newf2, label="f2trans")
                newf2 += self.f2func[symbol](newE)
                if self.debug:
                    pl.plot(E, f2, label="f2old")
                    pl.plot(newE, newf2, label="f2+f2tab")
                #newE, newf2 = mskk.imag(E, thisf1)
                #        self.linearf1[symbol] = lambda x: f1[ileft] + (emax - emin)/(f1[iright] - f1[ileft]) * (x - emin)
                f2func = interpolate.interp1d(newE, newf2)
                newf2 = f2func(E[ind])
                #print newf2[0], f2[il], E[ind][0], E[il]
                #newf2 -= newf2[0] + (newf2[-1] - newf2[0]) / (E[ind][-1] - E[ind][0]) * (E[ind] - E[ind][0])
                #newf2 +=    f2[il] + (   f2[ir] -    f2[il]) / (E[ir] - E[il]) * (E[ind] - E[il])
                f2[ind] = newf2
            else:
                KK[symbol] = "f2"
                if not ay is None:
                    ay = ay.real - self.Z[symbol]
                """
                ancind = abs(E - np.array(ax)[:,np.newaxis]).argmin(-1)
                ay = f1[ancind]
                ax = E[ancind]
                """
                thisf2 = f2 - self.f2func[symbol](E)
                newE, newf1 = mskk.real(E, thisf2, ax, ay, corr=True)
                if self.debug:
                    import pylab as pl
                    pl.plot(E, f2, label="f2")
                    pl.plot(E, self.f2func[symbol](E), "--", label="f2tab")
                    pl.plot(E, thisf2, label="f2 - f1tab")
                    pl.plot(newE, newf1, label="f1trans")
                newf1 += self.f1func[symbol](newE)
                if self.debug:
                    pl.plot(E, f1, label="f1old")
                    pl.plot(newE, newf1, label="f1+f1tab")
                #newE, newf1 = mskk.real(E, f2)
                f1func = interpolate.interp1d(newE, newf1)
                newf1 = f1func(E[ind])
                #print newf1[0], f1[il], E[ind][0], E[il]
                #newf1 -= newf1[0] + (newf1[-1] - newf1[0]) / (E[ind][-1] - E[ind][0]) * (E[ind] - E[ind][0])
                #newf1 +=    f1[il] + (   f1[ir] -    f1[il]) / (E[ir] - E[il]) * (E[ind] - E[il])
                f1[ind] = newf1 # f1func(E[ind])
            if self.debug:
                import pylab as pl
                pl.legend()
                pl.xlabel("energy / eV")
                pl.ylabel("f / electrons")
                #pl.ylim(ymax=90, ymin=0)
                pl.grid()
                pl.show()
            self.f[symbol.name] = f1 + 1j * f2
        self.Isim = abs(self.Ffunc.dictcall(self.f))**2
        err = ((self.DAFS[allind] - self.Isim[allind])**2).sum()
        return err, KK


