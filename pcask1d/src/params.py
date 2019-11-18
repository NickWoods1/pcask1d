# Distributed under the terms of the MIT License

"""
This module implements the parameters class, which carries the input model for the computation.
An object of type Parameters is designed to be immutable and lead to a unique output.
"""

import numpy as np


class Parameters:
    """
    The use of this class is to generate a unique input for a DFT calculation.
    This includes parameters set by the user (e.g. scf_tol), and any derived quantities
    (e.g. num_planewaves).
    """

    def __init__(self, **kwargs):

        # Level of approximation used
        self._method = kwargs.get('method', 'dft')

        # Size of real space cell
        self._cell = kwargs.get('cell', 20)
        self._num_planewaves = kwargs.get('num_planewaves', 1001)
        self._k_point_spacing = kwargs.get('kpoint_spacing', 0.2)

        # List of species + position
        self._species = kwargs.get('species', ['Li'])
        self._positions = kwargs.get('positions', [0])

        # Define each element with corresponding atomic number
        self._element_charges = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5,
                                 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10}

        # SCF parameters
        self._scf_tol = kwargs.get('scf_tol', 1e-10)
        self._scf_history_length = kwargs.get('scf_history_length', 10)
        self._scf_step_length = kwargs.get('scf_step_length', 1)
        self._scf_temperature = kwargs.get('scf_temperature', 300)

        # Pseudopotential
        self._soft = 0.1

        # Have v_ext specified by atoms or given explicitly
        self._manual_v_ext = kwargs.get('manual_v_ext', None)

        # Sanity checks
        assert len(self._species) == len(self._positions), 'Each element requires a unique position.'
        assert abs(max(self._positions)) <= self._cell, 'All elements must lie within the primitive unit cell.'
        assert self._num_planewaves % 2 != 0, 'Number of plane-waves must be odd.'

        if self._method not in ['h', 'hf', 'dft']:
            raise RuntimeError('Chosen method of {} is not implemented'.format(self._method))

    def smearing_scheme(self, energy):
        """ Ansatz for smearing the occupancies to prevent
        occupancy-induced instability in SCF iterations """
        return 1 / (np.exp(energy / self._scf_temperature) + 1)

    @property
    def element_charges(self):
        return self._element_charges

    @property
    def positions(self):
        return self._positions

    @property
    def species(self):
        return self._species

    @property
    def soft(self):
        r""" Defines the parameter c for the softened Coulomb potential :math:`\frac{1}{|x-x'| + c}` """
        return self._soft

    @property
    def num_planewaves(self):
        r""" Number of plane-waves :math:`\{ | G \rangle \}` in the basis set expansion """
        return self._num_planewaves

    @property
    def cell(self):
        r""" Size of the periodic unit cell, a, such that :math:`\Sigma \in [-a, a]` """
        return self._cell

    @property
    def scf_tol(self):
        r""" Convergence tolerance for fluctuations in the residual norm :math:`||R||` """
        return self._scf_tol

    @property
    def scf_history_length(self):
        r""" Number of iterations to include in the iterative SCF history of densities if
          using a quasi-Newton method for SCF convergence """
        return self._scf_history_length

    @property
    def scf_step_length(self):
        r""" Damping parameter :math:`\alpha \in [0,1)` applied to SCF steps """
        return self._scf_step_length

    @property
    def num_electrons(self):
        """ Number of electrons such that charge neutrality is enforced """
        num_electrons = 0
        for i in range(len(self._species)):
            num_electrons += self._element_charges[self._species[i]]
        return num_electrons

    @property
    def realspace_grid(self):
        """ Grid points in the delta function (real-space) basis set """
        return np.linspace(-self._cell, self._cell, self._num_planewaves)

    @property
    def big_realspace_grid(self):
        """ Real space grid with double the sampling """
        return np.linspace(-self._cell, self._cell, 2*self._num_planewaves)

    @property
    def planewave_grid(self):
        r""" Plane-wave frequences for plane-waves that fit in
         the unit cell: :math:`G = \frac{2 \pi n}{R}` """
        N = int((self._num_planewaves - 1) / 2)

        pw_frequencies_positive = [np.pi*n / self._cell
                          for n in
                          range(0, N+1)]

        pw_frequencies_negative = [np.pi*n / self._cell
                          for n in
                          range(-N, 0)]
        # Return frequencies ordered using same ordering as numpy's FFT
        return np.asarray(pw_frequencies_positive + pw_frequencies_negative)

    @property
    def big_planewave_grid(self):
        r""" Plane-wave frequences up to 2*G_max """
        N = self._num_planewaves

        pw_frequencies_positive = [np.pi*n / self._cell
                          for n in
                          range(0, N+1)]

        pw_frequencies_negative = [np.pi*n / self._cell
                          for n in
                          range(-N, 0)]

        return np.asarray(pw_frequencies_positive + pw_frequencies_negative)

    @property
    def v_ext(self):
        """ External potential: either an atomic potential (Coulomb) or a given functional form """
        if self._manual_v_ext is not None:
            return np.fft.fft(self._manual_v_ext(x=self.realspace_grid))
        else:
            v_ext = np.zeros(self._num_planewaves)
            num_atoms = len(self._species)
            # Add a Coulomb potential at each atomic position scaled with Z
            for i in range(num_atoms):
                charge = self._element_charges[self._species[i]]
                v_ext += self.coulomb(charge, self._positions[i])
            return np.fft.fft(v_ext)

    @property
    def k_points(self):
        r""" MP k-point grid within the first BZ: :math:`k \in [\frac{-\pi}{a}, \frac{\pi}{a}]` """
        return np.linspace(-np.pi / self._cell, np.pi / self._cell, 1 / self._k_point_spacing)

    def coulomb(self, charge: int, position: float) -> np.ndarray:
        """ The external potential of an ion with a given charge (int) and position (float)
         regularised about the core with a softening parameter """
        return -charge / (abs(self.realspace_grid - position) + self._soft)

