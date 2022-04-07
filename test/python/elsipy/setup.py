from setuptools import find_packages
from numpy.distutils.core import setup, Extension
import subprocess
import os
import re
import warnings

import numpy as np
import mpi4py

def runcommand(cmd):
    process = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, universal_newlines=True)
    c = process.communicate()

    if process.returncode != 0:
        raise Exception("Something went wrong whilst running the command: %s" % cmd)

    return c[0]

def whichmpi():
    # Figure out which MPI environment this is
    import re
    try:
        mpiv = runcommand('mpirun -V')

        if re.search('Intel', mpiv):
            return 'intelmpi'
        elif re.search('Open MPI', mpiv):
            return 'openmpi'
    except:
        return 'mpich'

    warnings.warn('Unknown MPI environment.')
    return None

def whichscalapack():
    # Figure out which Scalapack to use
    if 'MKLROOT' in os.environ:
        return 'intelmkl'
    else:
        return 'netlib'

mpiversion = whichmpi()

scalapackversion = whichscalapack()



## Find the MPI arguments required for building the modules.
if mpiversion == 'intelmpi':
    # Fetch command line, convert to a list, and remove the first item (the command).
    intelargs = runcommand('mpicc -show').split()[1:]
    mpilinkargs = intelargs
    mpicompileargs = intelargs
elif mpiversion == 'openmpi':
    # Fetch the arguments for linking and compiling.
    mpilinkargs = runcommand('mpicc -showme:link').split()
    mpicompileargs = runcommand('mpicc -showme:compile').split()
elif mpiversion == 'mpich':
    link_info = runcommand('mpicc -link_info')
    #mpilinkargs = re.sub(r'^\w+\s', '', link_info).split()
    mpilinkargs = link_info.split()[1:]
    compile_info = runcommand('mpicc -compile_info')
    #mpicompileargs = re.sub(r'^\w+\s', '', compile_info).split()
    mpicompileargs = compile_info.split()[1:]
else:
    raise Exception("MPI library unsupported. Please modify setup.py manually.")

## Find the Scalapack library arguments required for building the modules.
if scalapackversion == 'intelmkl':
    # Set library includes (taking into account which MPI library we are using)."
    scl_lib = ['mkl_scalapack_lp64', 'mkl_rt', 'mkl_blacs_'+mpiversion+'_lp64', 'mkl_intel_lp64', 'mkl_sequential', 'mkl_core',      'iomp5', 'pthread']
    scl_libdir = [os.environ['MKLROOT']+'/lib/intel64' if 'MKLROOT' in os.environ else '']
elif scalapackversion == 'netlib':
    scl_lib = ['scalapack', 'gfortran']
    scl_libdir = [ os.path.dirname(runcommand('gfortran -print-file-name=libgfortran.a')) ]
else:
    raise Exception("Scalapack distribution unsupported. Please modify setup.py manually.")

elsi_libdir=["/home/yy244/elsi/elsipy/elsi_interface/build_2/lib"]
elsi_lib = ["elsi", "elpa", "fortjson", "MatrixSwitch", "NTPoly", "OMM"]

## Try and decide whether to use Cython to compile the source or not.
try:
    from Cython.Build import cythonize
    HAVE_CYTHON = True
except ImportError as e:
    warnings.warn("Cython not installed.")
    HAVE_CYTHON = False

blacs_ext = Extension('elsipy.blacs', ["elsipy/blacs.pyx"],
                      include_dirs=['.',  np.get_include(), mpi4py.get_include()],
                      library_dirs=scl_libdir, libraries=scl_lib,
                      extra_compile_args=mpicompileargs,
                      extra_link_args=mpilinkargs,
                      define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])

elsi_ext = Extension('elsipy.elsi', ["elsipy/elsi.pyx"],
                      include_dirs=['.',"/home/yy244/elsi/elsipy/elsi_interface/src/include",  np.get_include(), mpi4py.get_include()],
                      library_dirs=elsi_libdir, libraries=elsi_lib,
                      extra_compile_args=mpicompileargs,
                      extra_link_args=mpilinkargs,
                      define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])

exts=[blacs_ext, elsi_ext]

#exts = cythonize(exts, include_path=['.', np.get_include(), mpi4py.get_include()],compiler_directives={'language_level' : "3"})
exts = cythonize(exts, include_path=['.', np.get_include(), mpi4py.get_include()])

setup(
    name='elsipy',
    author='Yi Yao',
    description='Python bindings for ELSI.',
    url='http://github.com/',
    license='BSD',
    packages=find_packages(),
    ext_modules=exts
)
