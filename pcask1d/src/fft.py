# Distributed under the terms of the MIT License.

"""
Implements discrete Fourier transform over period params.cell \in [-a,a]
and frequencies defined by G = n*pi/a
"""

import numpy as np


class Fourier:
    """
    Interface for switching between spaces for a given params
    """
    def __init__(self, params):
        self._planewave_grid = params.planewave_grid
        self._realspace_grid = params.realspace_grid

    @staticmethod
    def fft(function_real):
        """ Fourier transform: f(G) = (1/N) \sum f(x) e^{2pi i G x}"""
        return np.fft.fft(function_real)

    @staticmethod
    def ifft(function_fourier):
        """ Inverse Fourier transform """
        return np.fft.ifft(function_fourier)
