# Distributed under the terms of the MIT License.

"""
Module that defines a single-particle wavefunction and its various properties
"""

import numpy as np
from .density import Density


class Wavefunction:
    """ This class represents a single-particle wavefunction :math:`\phi_{ik}^\sigma (G)`
    with associated single-particle energy :math:`\varepsilon` and occupancy :math:`f_i`.
    Each wavefunction has a spin :math:`\sigma`, k-point :math:`k`, and band index :math:`i`"""

    def __init__(self, params, **kwargs):

        # Number of plane-waves
        N = len(params.planewave_grid)

        self._pw_coefficients = kwargs.get('pw_coefficients', np.zeros(N, dtype=complex))
        self._energy = kwargs.get('energy', 0)
        self._k_point = kwargs.get('k_point', 0)
        self._spin = kwargs.get('spin', 0)
        self._band_index = kwargs.get('band_index', 0)
        self._occupancy = kwargs.get('occupancy', 1)

    def __str__(self):
        return 'Band {0} of system Hamiltonian with k-point {1}'.format(self._band_index, self._k_point)

    @property
    def pw_coefficients(self):
        """ The plane-wave coefficients of the wavefunction """
        return self._pw_coefficients

    @pw_coefficients.setter
    def pw_coefficients(self, coeffs: np.ndarray) -> np.ndarray:
        self._pw_coefficients = coeffs

    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, energy: float):
        self._energy = energy

    @property
    def k_point(self):
        return self._k_point

    @k_point.setter
    def k_point(self, k_point: float):
        self._k_point = k_point

    @property
    def band_index(self):
        return self._band_index

    @band_index.setter
    def band_index(self, band_index: int):
        self._band_index = band_index

    @property
    def occupancy(self):
        return self._occupancy

    @occupancy.setter
    def occupancy(self, occ: float):
        self._occupancy = occ

    def get_density(self):
        """ Obtain the (occupancy weighted) single-particle density
         corresponding to a single-particle wavefunction """
        #TODO: put wavefunction on big grid, pad with zeros, then FFT.
        wavefunction = np.fft.ifft(self._pw_coefficients)
        return Density(coeffs=self._occupancy * wavefunction.conj() * wavefunction)
