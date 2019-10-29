import argparse
import numpy as np
from pcask1d.src.input import parameters

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

        # Construct parameters object
        params = parameters()

        kpoints = params.get_kpoints()
        planewave_frequencies = params.get_planewave_frequencies()

        hamiltonian = ham()

        wavefunction = []
        for k in kpoints:

            ham = hamiltonian(k)

            # Construct Hamiltonian
            # wavefunction.append(hamiltonian.get_eigenfunctions())
            # wavefunction[i].set_kpoint(k)
            # wavefunction[i].set_
            #
            # Diag Hamiltonian

            # Store orbitals

        # Construct density (over all kpoints)

        # do density mixing



        # Find g.s. wavefunctions, energy, and density
        #wavefunctions, total_energy, density = minimise_energy_hf(params)

        print('hello')
