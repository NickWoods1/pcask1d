# Distributed under the terms of the MIT License.:w

"""
The file defines a particle density (probability distribution) and its various features/uses.
"""

import numpy as np


class Density:
    """ The class representing a particle density """

    def __init__(self):

        self.coefficients = np.zeros(10)

        assert min(self.coefficients) >= 0, "Negative density region exists, aborting."

    def __add__(self, density):
        return self.coefficients + density.coefficients


    # Check norm of the density
    def norm(self):
        return None