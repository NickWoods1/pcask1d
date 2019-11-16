# Distributed under the terms of the MIT License.

"""
This module contains the Hamiltonian. Each Hamiltonian
is defined for a given k-point and density, with methods to
extract various quantities such as total energy, wavefunctions, and so on.
"""

import numpy as np
import scipy as sp
from scipy.sparse.linalg import eigs
import warnings


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
                                     + np.diag(self.v_h(params)) \
                                     + np.diag(self.v_xc(params))
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
            return np.linalg.eigh(self.representation(params))
        else:
            return eigs(self.representation(params), num_states, which='SM')

    def kinetic(self, params):
        r""" Kinetic operator in Fourier space: :math:`\hat{T} = \frac{1}{2} |G+k|^2`

        Parameters:
            * params (Parameters): input model for the system
         """

        return 0.5*abs(params.planewave_grid + self._k_point)**2

    def v_xc(self, params):
        """ Exchange-correlation potential

        Parameters:
            * params (Parameters): input model for the system
        """

        N = len(params.planewave_grid)
        return np.zeros(N)

    def v_h(self, params):
        r""" The Hartree potential: :math:`v_h = \frac{4 \pi \rho(G)}{|x-x'|}`

        Parameters:
            * params (Parameters): input model for the system
        """

        G = params.planewave_grid
        G[0] = 1
        return 4.0*np.pi*self._density / G**2

    def v_ext(self, params):
        r""" The non-local external potential operator:
        :math:`v_{ext}(G-G')` from :math:`\mathcal{F}(v_{ext}(x))`

        Parameters:
            * params (Parameters): input model for the system
        """

        N = len(params.planewave_grid)
        v_ext_operator = np.zeros((N,N), dtype=complex)
        v_ext = params.v_ext
        for i in range(N):
            j = 0
            while j <= i:
                v_ext_operator[i,j] = v_ext[i - j]
                j += 1
        v_ext_operator += v_ext_operator.T.conj() #FIX DIAGONAL
        return v_ext_operator

    def total_energy(self, wavefunctions=None):
        r""" Total energy of a Hamiltonian: :math:`E = \langle \Psi | H | \Psi \rangle`"""
        pass
