PEXSI: Pole EXpansion and Selected Inversion
============================================

The Pole EXpansion and Selected Inversion (PEXSI) method is a fast method for
electronic structure calculation based on Kohn-Sham density functional theory.
It efficiently evaluates certain selected elements of matrix functions, e.g.,
the Fermi-Dirac function of the KS Hamiltonian, which yields a density matrix.
It can be used as an alternative to diagonalization methods for obtaining the
density, energy and forces in electronic structure calculations. The PEXSI
library is written in C++, and uses message passing interface (MPI) to
parallelize the computation on distributed memory computing systems and achieve
scalability on more than 10,000 processors.

From numerical linear algebra perspective, the PEXSI library can be used as a
general tool for evaluating certain selected elements of a matrix function, and
therefore has application beyond electronic structure calculation as well.

The documentation of PEXSI is compiled by Sphinx hosted on

https://pexsi.readthedocs.io

For installation instructions please (be patient for dependecies) see

https://pexsi.readthedocs.io/en/latest/install.html
