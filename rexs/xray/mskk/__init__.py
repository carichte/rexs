from core import *
from digest import *

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

