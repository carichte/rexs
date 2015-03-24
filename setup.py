#!/usr/bin/env python
from numpy.distutils.core import setup, Extension
import numpy
import sys

if len(sys.argv)<2:
    print("see install.txt for installation instructions.")


setup( name = "rexs", 
       version = "0.2",
       ext_modules = [Extension("rexs.xray.deltaf.deltaf", ["rexs/xray/deltaf/deltaf.f"]), 
                      Extension("rexs.tools.rebin.librebin", ["rexs/tools/rebin/rebin.c"])],
       packages = ["rexs", 
                   "rexs.xray", 
                   "rexs.xray.mskk", 
                   "rexs.xray.deltaf", 
                   "rexs.io", 
                   "rexs.tools",
                   "rexs.tools.rebin"
                   ],
       package_data={"rexs": ["xray/deltaf/*.npz",
                              "xray/elementdata.sqlite"]
                     },
       author = "Carsten Richter", 
       author_email = "carsten.richter@desy.de",
       description = "rexs - toolkit for evaluation of resonant x-ray "
                     "scattering measurements.",
       long_description = """
                             rexs - toolkit for evaluation of resonant x-ray 
                             scattering measurements
                          """
     )

