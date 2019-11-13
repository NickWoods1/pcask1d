"""
Single-particle wavefunctions and their various properties
"""

import numpy as np


class Wavefunction:

    def __init__(self, params, **kwargs):

        # Number of plane-waves
        N = len(params.planewave_grid)

        self.pw_coefficients = kwargs.get('pw_coefficients', np.zeros(N, dtype=complex))
        self.energy = kwargs.get('energy', 0)
        self.k_point = kwargs.get('k_point', 0)
        self.spin = kwargs.get('spin', 0)
        self.band_index = kwargs.get('band_index', 0)

    def __str__(self):
        return 'Band {0} of system Hamiltonian with k-point {1}'.format(self.band_index, self.k_point)
