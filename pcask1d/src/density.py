# Distributed under the terms of the MIT License.

"""
The module defines a particle density (probability distribution)
and its various features/uses.
"""

import numpy as np


class Density:
    """ The class representing a particle density """

    def __init__(self, params, **kwargs):

        self._coefficients = kwargs.get('coeffs', self.initial_guess(params))
        assert min(self._coefficients) >= 0, "Negative density region exists, aborting."

    def __add__(self, density):
        return self._coefficients + density.coefficients

    @property
    def coefficients(self):
        return self._coefficients

    @coefficients.setter
    def coefficients(self, coeffs: np.ndarray):
        self._coefficients = coeffs

    @staticmethod
    def initial_guess(params):
        """ Initial guess for the density as overlapping Gaussians of charge """
        density = np.zeros(2*params.num_planewaves)
        for i in range(len(params.species)):
            charge = params.element_charges[params.species[i]]
            density += charge*np.exp(-(params.big_realspace_grid - params.positions[i])**2)

        density = np.fft.fft(density)
        density[0] = params.num_electrons
        return density

    def norm(self):
        r""" The L1 norm of the density:

        .. math::

            \int_{\Sigma} \rho(x) dx = N
        """
        return self._coefficients[0]
