import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


class Hamiltonian:

    # Default constructor
    def __init__(self, density, **kwargs):

        # Default gamma point
        self.k_point = kwargs.pop('k_point',0)
        self.density = density

    def representation(self, params):
        """ Get the representation (coefficients) of the Hamiltonian -- plane-wave basis"""
        hamiltonian_representation = self.kinetic(params) \
                                     + self.v_xc(params, self.density) \
                                     + self.v_ext(params) \
                                     + self.v_h(params, self.density)
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
        return 1 / abs(params.planewave_grid + self.k_point)

    def v_xc(self, params, density):
        return np.zeros(params.num_planewaves)

    def v_h(self, params, density):
        return np.zeros(params.num_planewaves)

    def v_ext(self, params):
        # REQUIRES non-local "projectors" for the FFT?
        return params.v_ext
