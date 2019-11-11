"""
Implements discrete Fourier transform over period params.cell \in [-a,a]
and frequencies defined by G = n*pi/a
"""

import matplotlib.pyplot as plt
import numpy as np


class Fourier:
    """
    Interface for switching between spaces for a given params
    """
    def __init__(self, params):
        self._planewave_grid = params.planewave_grid
        self._realspace_grid = params.realspace_grid

    @staticmethod
    def fft(params, function_real):
        """ Discrete Fourier transform: f(G) = (1/N) \sum f(x) e^{2pi i G x}"""

        N = len(params.realspace_grid)
        function_fourier = np.zeros(N, dtype=complex)
        for G in range(N):
            function_fourier[G] = (1/(N))*np.dot(function_real,
                                                 np.exp(1.0j*params.planewave_grid[G]*params.realspace_grid))

        return function_fourier

    @staticmethod
    def ifft(params, function_fourier):
        """ Inverse discrete Fourier transform: f(G) = \sum f(x) e^{-2pi i G x}"""

        N = len(params.realspace_grid)
        function_real = np.zeros(N)
        for x in range(N):
                function_real[x] = np.dot(function_fourier,
                                          np.exp(-1.0j*params.realspace_grid[x]*params.planewave_grid))

        return function_real
