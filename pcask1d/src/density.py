"""
Implements class that defines a particle density (probability distribution)
"""

import numpy as np


class Density:

    def __init__(self):

        self.coefficients = np.zeros(10)

        assert min(self.coefficients) >= 0, "Negative density region exists, aborting."

    def __add__(self, density):
        return self.coefficients + density.coefficients


    # Check norm of the density
    def norm(self):
        return None