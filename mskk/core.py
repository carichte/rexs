#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Mo 7. Mar 16:44:56 CEST 2013
# Computer: haso227r 
# System: Linux 3.13.0-30-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
import numpy as np
from scipy import ndimage, integrate

__doc__ =  """
    The Module provides functions for performing a kramers-kronig transform
    on given data of the real or imaginary part of the susceptibility or 
    similar optical properties.
    It supports multiple substractive kramers-kronig transforms according 
    to the work of e.g. Lucarini et al. (J. Chem. Phys. 119, 11095 (2003))
    Therefore one has to provide anchor points of the quantity which is to
    be calculated.
    
    Functions:
        
        real : get the real part from given data of the imaginary part
        imag : get the imaginary part from given data of the real part
"""



def real(om, imag, anc_om=None, anc_re=None, verbose=True, corr=True):
    """
        Calculate the real part of the susceptibility (or similar) out of
        the imaginary part using the kramers-kronig transform and optionally
        a set of known anchor points.
        
        Further reading: J. Chem. Phys. 119, 11095 (2003)
                         http://dx.doi.org/10.1063/1.1623477
        This is the special case where alpha=0 and n=1. 
        (No moments and only linear optics)
    """
    if anc_om==None or anc_re==None:
        Q=0
        anc_om=np.array(())
        anc_re=np.array(())
    else:
        om = np.array(om, dtype=float,ndmin=1)
        imag = np.array(imag, dtype=float,ndmin=1)
        anc_om = np.array(anc_om, dtype=float,ndmin=1)
        anc_re = -np.array(anc_re, dtype=float,ndmin=1)
        assert np.ndim(imag)==1, \
            "imag must be 1 dimensional"
        assert np.shape(anc_om)==np.shape(anc_re) and np.ndim(anc_re)==1, \
            "shape mismatch: anc_om and anc_re must have equal length and 1 dimension"
        Q = len(anc_om)
    
    
    om_prim = om.copy()[:,np.newaxis]
    dom = np.unique(np.diff(om))
    if len(dom)>1: 
        if (abs(dom/dom[0]-1)>1e-6).any():
            raise ValueError("Need equally spaced data.")
        else:
            dom = dom[0]
    om = np.append(om - dom/2., om[-1]+dom/2.)
    
    ind = abs(om - anc_om[:,np.newaxis]).argmin(1)
    if ((abs(anc_om - om[ind])/om[ind]) > 0.01).any() and verbose==True: 
        print("Warning: Anchor point far from measured data point")
    ind2 = (ind>0) * (ind<(len(om)-1)) + ((abs(anc_om - om[ind])/om[ind]) < 0.01)
    anc_om[ind2] = om[ind[ind2]]
    
    result = np.zeros(len(om))
    for j in range(Q):
        numerator = (om**2 - (np.delete(anc_om, j)**2)[:,np.newaxis]).prod(0)
        denominator = (anc_om[j]**2 - np.delete(anc_om, j)**2).prod()
        result += numerator/denominator * anc_re[j]
    
    # now the integral:
    imag = imag[:,np.newaxis]
    
    prefac = (om**2 - (anc_om**2)[:,np.newaxis]).prod(0) * 2/np.pi
    denominator1 = (om_prim**2 - anc_om**2).prod(1)[:,np.newaxis]
    denominator = ((om_prim**2 - om**2) * denominator1)
    integrand = om_prim * imag / denominator
    
    #result += prefac * integrate.simps(integrand, dx=dom, axis=0)
    """ 
    # Simpson manuell
    ind = np.arange((len(om_prim)-1)/2)*2
    I = dom/6*(integrand[ind] \
             + integrand[ind+2] \
           + 4*integrand[(ind+1)]).sum(0)
    result += prefac * I
    """
    #I=(h/6)*(funktion(a+i*h)+4*funktion(a+(i+0.5)*h)+funktion(a+(i+1)*h))
    
    # Simpson does not work
    result += prefac * integrate.trapz(integrand, dx=dom, axis=0)
    
    if len(anc_om)==0 and corr:
        #print "adding corr"
        omleft = om_prim[0].item() #- dom
        omrigh = om_prim[-1].item() #+ dom
        result -= prefac * (imag[0]  * np.arctanh((omleft/om).astype(complex)).real \
                          - imag[-1] * np.arctanh((omrigh/om).astype(complex)).real)
    return om, -result



def imag(om, real, anc_om=None, anc_im=None, verbose=True, corr=True):
    """
        Calculate the imaginary part of the susceptibility (or similar) out of
        the real part using the kramers-kronig transform and optionally
        a set of known anchor points.
        
        Further reading: J. Chem. Phys. 119, 11095 (2003)
                         http://dx.doi.org/10.1063/1.1623477
        This is the special case where alpha=0 and n=1. 
        (No moments and only linear optics so far)
    """
    if anc_om==None or anc_im==None:
        Q=0
        anc_om=np.array(())
        anc_im=np.array(())
    else:
        om = np.array(om, dtype=float,ndmin=1)
        real = np.array(real, dtype=float,ndmin=1)
        anc_om = np.array(anc_om, dtype=float,ndmin=1)
        anc_im = np.array(anc_im, dtype=float,ndmin=1)
        assert np.ndim(real)==1, \
            "real must be 1 dimensional"
        assert np.shape(anc_om)==np.shape(anc_im) and np.ndim(anc_im)==1, \
            "shape mismatch: anc_om and anc_im must have equal length and 1 dimension"
        Q = len(anc_om)
    
    
    om_prim = om.copy()[:,np.newaxis]
    dom = np.unique(np.diff(om))
    
    if len(dom)>1: 
        if (abs(dom/dom[0]-1)>1e-6).any():
            raise ValueError("Need equally spaced data.")
        else:
            dom = dom[0]
    om = np.append(om - dom/2., om[-1]+dom/2.)
    
    ind = abs(om - anc_om[:,np.newaxis]).argmin(1)
    if ((abs(anc_om - om[ind])/om[ind]) > 0.01).any() and verbose==True: 
        print("Warning: Anchor point far from measured data point")
    ind2 = (ind>0) * (ind<(len(om)-1)) + ((abs(anc_om - om[ind])/om[ind]) < 0.01)
    anc_om[ind2] = om[ind[ind2]]
    
    result = np.zeros(len(om))
    for j in range(Q):
        numerator = (om**2 - (np.delete(anc_om, j)**2)[:,np.newaxis]).prod(0)
        denominator = (anc_om[j]**2 - np.delete(anc_om, j)**2).prod()
        result += numerator/denominator * anc_im[j] / anc_om[j]
    
    # now the integral:
    real = -real[:,np.newaxis]
    
    prefac = (om**2 - (anc_om**2)[:,np.newaxis]).prod(0) * 2/np.pi
    denominator1 = (om_prim**2 - anc_om**2).prod(1)[:,np.newaxis]
    denominator = ((om_prim**2 - om**2) * denominator1)
    integrand = real / denominator
    
    # Simpson does not work
    result -= prefac * integrate.trapz(integrand, dx=dom, axis=0)
    result *= om
    if len(anc_om)==0 and corr:
        omleft = om_prim[0].item() #- dom
        omrigh = om_prim[-1].item() #+ dom
        result += prefac * (real[0]  * np.arctanh((omleft/om).astype(complex)).real \
                          - real[-1] * np.arctanh((omrigh/om).astype(complex)).real)
    return om, result

