import os
import sys
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class CustomBuildExt(build_ext):
    def build_extension(self, ext):
        if ext.name == "rexs.xray.deltaf.deltaf":
            # The name you'll give to the compiled extension module
            so_name = "deltaf"
            
            # Absolute path to your Fortran source file:
            # e.g., mypkg/rexs/xray/deltaf/deltaf.f
            this_dir = os.path.abspath(os.path.dirname(__file__))
            deltaf_source = os.path.join(
                this_dir, "rexs", "xray", "deltaf", "deltaf.f"
            )
            
            # Create the build_temp directory if it doesn't exist
            build_temp = os.path.abspath(self.build_temp)
            os.makedirs(build_temp, exist_ok=True)

            # f2py command, referencing the absolute path to deltaf.f
            cmd = [
                sys.executable, "-m", "numpy.f2py", "-c", "-m", so_name,
                deltaf_source,         # absolute path here
                #"-O3",                 # or other compiler flags
            ]

            # Invoke f2py in the build directory:
            subprocess.check_call(cmd, cwd=build_temp)

            # f2py should produce a .so (or .pyd) file named deltaf.* there
            built_so = os.path.join(build_temp, f"{so_name}.so")
            if not os.path.exists(built_so):
                # On Windows, f2py might generate .pyd or .dll 
                # so you might need to search for the correct file extension:
                # e.g. 'deltaf.cp310-win_amd64.pyd'
                # or adapt logic as needed.
                for f in os.listdir(build_temp):
                    if f.startswith(so_name) and (f.endswith(".so") or f.endswith(".pyd") or f.endswith(".dll")):
                        built_so = os.path.join(build_temp, f)
                        break

            # Move the compiled extension into the final location
            ext_path = self.get_ext_fullpath(ext.name)
            self.copy_file(built_so, ext_path)

        else:
            # For your .c extension, use the normal build_ext
            super().build_extension(ext)

ext_modules = [
    # We'll build this one ourselves with f2py:
    Extension("rexs.xray.deltaf.deltaf", sources=[]),
    # Plain C extension:
    Extension("rexs.tools.rebin.librebin", sources=["rexs/tools/rebin/rebin.c"]),
]

setup(
    name="rexs",
    version="0.2",
    author = "Carsten Richter", 
    author_email = "carsten.richter@desy.de",
    description = "rexs - toolkit for evaluation of resonant x-ray "
                     "scattering measurements.",
    cmdclass={"build_ext": CustomBuildExt},
    ext_modules=ext_modules,
    packages=[
        "rexs",
        "rexs.xray",
        "rexs.xray.mskk",
        "rexs.xray.deltaf",
        "rexs.io",
        "rexs.tools",
        "rexs.tools.rebin",
    ],
    package_data={"rexs": [
        "xray/deltaf/*.npz",
        "xray/elementdata.sqlite",
    ]},
    install_requires=[
        "numpy",
        "scipy",
        "meson",
        "ninja",
        "six"
        ],
)
