# Distributed under the terms of the MIT License.

"""
The module defines a particle density (probability distribution)
and its various features/uses.
"""

import numpy as np


class Density:
    """ The class representing a particle density """

    def __init__(self):

        self.coefficients = np.zeros(10)

        assert min(self.coefficients) >= 0, "Negative density region exists, aborting."

    def __add__(self, density):
        return self.coefficients + density.coefficients

    def norm(self):
        """ The L1 norm of the density
        ..math::

            \int_{\Sigma} \rho(x) dx = N
        """
        pass
