*******
### pCASK1d
*******

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



`pcask1d` is a python package designed to compute the ground
state energy and density of a molecular system in one-dimension using
various approximations. `pcask1d` is periodic
 and uses a plane-wave basis set. 

*****************
#### Functionality
*****************

Given an external potential, specified manually (e.g. harmonic well), or using the in-built utility
to specify atomic species & positions:

`dft` generates the Kohn-Sham self-consistent single-particle
wavefunctions, total energy, and density, given an exchange-correlation potential (currently, LDA).

`h` generates Hartree self-consistent single-particle 
wavefunctions, total energy, and density.
