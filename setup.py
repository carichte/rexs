#!/usr/bin/env python
from distutils.core import setup, Extension
import numpy
import sys

if len(sys.argv)<2:
    print("see install.txt for installation instructions.")

setup( name = "pyrexs", 
       version = "0.1",
       packages = ["mskk", "evaluationtools", "evaluationtools", "materials"],
       package_data={'materials': ['cif/*']},
       author = "Carsten Richter", 
       author_email = "carsten.richter@desy.de",
       description = "pyrexs - toolkit for evaluation of resonant x-ray scattering measurements.",
       long_description = """
                             pyrexs - toolkit for evaluation of resonant x-ray scattering measurements
                          """
     )

