elsipy provide a general parallel python interface to ELSI via mpi4py.

Steps to install elsipy, I use miniconda as package manager

1) conda create -n elsi_env python=3.9
2) conda activate elsi_env
3) pip install numpy scipy cython nose
4) conda install mpich scalapack mpi4py
5) After modify the paths in setup.py; python setup.py install

An example is given in test.py
