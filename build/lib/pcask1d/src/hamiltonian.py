"""
Class that implements the Hamiltonian object. Defined for a given k-point and density, with methods to
extract various quantities such as total energy, wavefunctions, and so on.
"""

import numpy as np
import scipy as sp


class Hamiltonian:

    # Default constructor
    def __init__(self, density, **kwargs):

        # For which k-point is H constructed? Default \gamma point.
        self._k_point = kwargs.pop('k_point',0)

        # Density associated with the Kohn-Sham Hamiltonian
        self._density = density

    def representation(self, params):
        """ Get the representation (coefficients) of the Hamiltonian -- plane-wave basis"""
        hamiltonian_representation = self.kinetic(params) \
                                     + self.v_ext(params) \
                                     + self.v_h(params) \
                                     + self.v_xc(params)
        return hamiltonian_representation

    def eigendecomposition(self, params, num_states='all'):
        """ Get the num_states lowest lying eigenvectors and eigenvalues of the Hamiltonian """
        if num_states == 'all':
            return np.linalg.eigh(self.representation(params))
        else:
            return sp.linalg.eigsh(self.representation(params), num_states)

    def total_energy(self, wavefunctions=None):
        """ Total energy of a Hamiltonian. (Optional: wavefunctions to avoid recalculating). """
        # if wvfns not None
        # get energy as sum over wvfn eigenenergies
        # else
        # get energy by diag then energy.

    def kinetic(self, params):
        """ Kinetic operator in Fourier space, 1/|G+k|^2"""
        return 0.5*abs(params.planewave_grid + self._k_point)**2

    def v_xc(self, params):
        return np.zeros(params.num_planewaves)

    def v_h(self, params):
        """ Hartree potential FT(rho / |x-x'|)"""
        G = params.planewave_grid
        G[0] = 1
        return 4.0*np.pi*self._density / G**2

    def v_ext(self, params):
        return params.v_ext
