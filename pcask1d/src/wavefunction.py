"""
Single-particle wavefunctions and their various properties
"""

import numpy as np


class Wavefunction:
    def __init__(self):

        self.pw_coefficients = 0
        self.kpoint = 0
        self.spin = 0
        self.band_index = 0
        self.energy = 0
