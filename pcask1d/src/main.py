import argparse
import numpy as np
import matplotlib.pyplot as plt
from pcask1d.src.params import Parameters
from pcask1d.src.hamiltonian import Hamiltonian

"""
Entry point for the requested action
"""

def main():

    __version__ = 0.1

    # Parser class for density2potential
    parser = argparse.ArgumentParser(
        prog='pcask1d',
        description='Find the ground state energy of  '
                    ' a periodic molecular system given an approximation',
        epilog='written by Nick Woods')

    # Specify arguments that the package can take
    parser.add_argument('--version', action='version', version='This is version {0} of cask1d.'.format(__version__))
    parser.add_argument('task', help='what do you want pcask1d to do: dft')

    args = parser.parse_args()

    # Code header
    print('    ------------------------')
    print('             pcask1d        ')
    print('    ------------------------')
    print('           Written by')
    print('           Nick Woods')
    print(' ')

    # Find the Kohn-Sham potential that generates a given reference density
    if args.task == 'get-groundstate':

        def v_ext(x):
            return 0.5*0.25**2*x**2

        # Construct parameters object
        params = Parameters(method='h',
                            species=['Li'],
                            positions=[0],
                            manual_v_ext=v_ext
                            )

        print(params.v_ext)

        #kpoints = params.get_kpoints()
        #planewave_frequencies = params.get_planewave_frequencies()

        # Construct density object with a given initial guess (part of the constructor)
        #density = particle_density(params)

        # Finds output density given an input density
        # wavefunction = []
        #for k in kpoints:

            # Create the Hamiltonian at the given k-point
            #H = hamiltonian(params, kpoint=k, density=density)

            # Returns a container of type wavefunction. Each with occupation, band index, etc.
            #wavefunctions = H.get_eigenfunctions(density)

            # Get density at k-point
            #for wavefunction in wavefunctions:
                #density += wavefunction.get_density()

        # Integrate over k-space to get the output density

        # Above, we just need F[rhoin] = rhoout, rest can be done in dm.py.
        # All done in Fourier space, with FFT's that transform between the two spaces (set this up
        # for periodic functions). Provide v_ext, initial guess density, in real space, represent in
        # periodic PW basis when needed. Need scheme to integrate over k-space when needed. (some interpolation?)

        # Can distribute over k-points if needed...

        # Iterative diag ---> explicit dense construction and dense diagonalisation
        # Compute all eigenstates (not just occupied) or just occupied (Lanszcos)
        # (Don't need to switch between spaces to apply operator)
        # Construct density in real space most efficient? (no?)
        # No pseudopotentials, all electron.
        # Use @property not get/set/del?
        # underscore to privatise variables that later can use get/set on

        # Guess density
        # While error > tol
        # Create H
        # Get output density
        # Mix densities