import argparse
import numpy as np
from p-cask1d.src.input import parameters

"""
Entry point for the requested action
"""

def main():

    __version__ = 0.1

    # Parser class for density2potential
    parser = argparse.ArgumentParser(
        prog='p-cask1d',
        description='Find the ground state energy of  '
                    ' a periodic molecular system given an approximation',
        epilog='written by Nick Woods')

    # Specify arguments that the package can take
    parser.add_argument('--version', action='version', version='This is version {0} of cask1d.'.format(__version__))
    parser.add_argument('task', help='what do you want cask1d to do: dft')

    args = parser.parse_args()

    # Code header
    print('    ------------------------')
    print('            p-cask1d        ')
    print('    ------------------------')
    print('           Written by')
    print('           Nick Woods')
    print(' ')

    # Find the Kohn-Sham potential that generates a given reference density
    if args.task == 'get-groundstate':

        # Construct parameters class
        params = parameters()

        # Find g.s. wavefunctions, energy, and density
        wavefunctions, total_energy, density = minimise_energy_hf(params)

