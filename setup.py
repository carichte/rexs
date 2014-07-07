#!/usr/bin/env python
from numpy.distutils.core import setup, Extension
import numpy
import sys

if len(sys.argv)<2:
    print("see install.txt for installation instructions.")

setup( name = "rexs", 
       version = "0.1",
       ext_modules = [Extension("deltaf.deltaf", ["deltaf/deltaf.f"])],
       packages = ["mskk", "evaluationtools", "materials", "deltaf"],
       package_data={'materials': ['cif/*'],
                        "deltaf": ["*.npz"]},
       author = "Carsten Richter", 
       author_email = "carsten.richter@desy.de",
       description = "rexs - toolkit for evaluation of resonant x-ray scattering measurements.",
       long_description = """
                             rexs - toolkit for evaluation of resonant x-ray scattering measurements
                          """
     )

