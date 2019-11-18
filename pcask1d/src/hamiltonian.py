# Distributed under the terms of the MIT License.

"""
This module contains the Hamiltonian. Each Hamiltonian
is defined for a given k-point and density, with methods to
extract various quantities such as total energy, wavefunctions, and so on.
"""

import numpy as np
import scipy as sp
import warnings
from scipy.sparse.linalg import eigs
from .wavefunction import Wavefunction


class Hamiltonian:
    """ Class that implements the Hamiltonian object """

    def __init__(self, density, **kwargs):
        r"""
        Defines (abstractly) the Hamiltonian operator for a given density and k-point from
        which all other quantities are derived.

        Parameters:
            * density (Density): density from which Hartree and XC potentials constructed

        Keyword arguments:
            * k-point: point on the reciprocal lattice for sampling the first BZ.
        """

        # For which k-point is H constructed? Default \gamma point.
        self._k_point = kwargs.get('k_point', 0)

        # Density associated with the Kohn-Sham Hamiltonian
        self._density = density

    def representation(self, params):
        r""" Construct the representation (coefficients) of the Hamiltonian
        in the plane-wave basis: :math:`\langle G | H[\rho] | G' \rangle`

        Parameters:
            * params (Parameters): input model for the system

        Output:
            * hamiltonian_representation (ndarray): the Hamiltonian matrix in plane-wave basis
        """

        hamiltonian_representation = np.diag(self.kinetic(params)) \
                                     + self.v_ext(params) \
                                     + self.v_h(params) \
                                     + self.v_xc(params)
        return hamiltonian_representation

    def eigendecomposition(self, params, num_states='all'):
        r""" Calculate the num_states lowest lying eigenvectors and eigenvalues of the Hamiltonian

        Parameters:
            * params (Parameters): input model for the system.
            * num_states (int): number of eigenvectors/eigenvalues to calculate.

        Output:
            * wavefunctions (Wavefunction): a container of type Wavefunction with
              the lowest num_states eigenvectors, band indices, etc..
        """

        if num_states > params.num_planewaves:
            num_states = 'all'
            warnings.warn('Requested num_states greater than Hamiltonian dimension -- calculating all eigenvectors.')

        if num_states == 'all':
            num_states = params.num_planewaves
            eigenvalues, eigenvectors = np.linalg.eigh(self.representation(params))
        else:
            eigenvalues, eigenvectors = eigs(self.representation(params), num_states, which='SM')

        wavefunctions = [Wavefunction for i in range(num_states)]
        for i, wavefunction in enumerate(wavefunctions):
            wavefunction.band_index = i
            wavefunction.pw_coefficients = eigenvectors[:,i]
            wavefunction.energy = eigenvalues[i]
            # TODO: Add support for partial occupancies here
            wavefunction.occupancy = 1 if i <= params.num_electrons else 0

        return wavefunctions

    def kinetic(self, params):
        r""" Kinetic operator in Fourier space: :math:`\hat{T} = \frac{1}{2} |G+k|^2`

        Parameters:
            * params (Parameters): input model for the system
         """

        N = params.num_planewaves
        return (0.5/N**2)*abs(params.planewave_grid + self._k_point)**2

    def v_xc(self, params):
        """ Exchange-correlation potential

        Parameters:
            * params (Parameters): input model for the system
        """

        return np.zeros((params.num_planewaves, params.num_planewaves))

    def v_h(self, params):
        r""" The Hartree potential in 1D

        Parameters:
            * params (Parameters): input model for the system
        """

        return self._density * np.fft.fft(1 / (abs(params.realspace_grid) + params.soft))

    def v_ext(self, params):
        r""" The non-local external potential operator:
        :math:`v_{ext}(G-G')` from :math:`\mathcal{F}(v_{ext}(x))`

        Parameters:
            * params (Parameters): input model for the system
        """

        return sp.linalg.toeplitz(params.v_ext)

    def total_energy(self, wavefunctions=None):
        r""" Total energy of a Hamiltonian: :math:`E = \langle \Psi | H | \Psi \rangle`"""
        pass
