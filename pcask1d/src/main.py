import argparse
import numpy as np
import matplotlib.pyplot as plt
from .params import Parameters
from .hamiltonian import Hamiltonian
from .wavefunction import Wavefunction
from .fft import Fourier

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
    if args.task == 'test':

        def v_ext(x):
            return 0.5*(0.25**2)*x**2#np.exp(-(1.5*x)**2)
            #return abs(x)#np.exp(-(1.5*x)**2)#0.5*(0.25**2)*x**2#np.exp(-(1.5*x)**2)

        # Construct parameters object
        params = Parameters(method='h',
                            species=['Li'],
                            positions=[0],
                            )

        density = v_ext(params.realspace_grid)
        hamiltonian = Hamiltonian(density)
        wavefunctions = [Wavefunction(params) for i in range(params.num_electrons)]

        eigenenergies, eigenfunctions = hamiltonian.eigendecomposition(params)

        for i, wavefunction in enumerate(wavefunctions):
            wavefunction.pw_coefficients = eigenfunctions[:,i]
            wavefunction.energy = eigenenergies[i]
            wavefunction.band_index = i
            wavefunction.k_point = 0


        #for wvfn in wavefunctions:
        #    wvfn.pw_coefficients = eigenfunction[:,i]

        #energies, wavefunctions = H.eigendecomposition(params)#, num_states=1)

        #plt.plot(wavefunctions[:,0])
        #plt.show()

      # Construct density object with a given initial guess (part of the constructor)
        #density = particle_density(params)
        # wvfn container: Finds output density given an input density
        # wavefunction = []
        #for k in kpoints:

            # Create the Hamiltonian at the given k-point
            #H = hamiltonian(params, kpoint=k, density=density)

            # Returns a container of type wavefunction. Each with occupation, band index, etc.
            #wavefunctions = H.get_eigenfunctions(density)

            # Get density at k-point
            #for wavefunction in wavefunctions:
                #density += wavefunction.get_density()

            #del H

