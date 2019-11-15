# Distributed under the terms of the MIT License

"""
Implements the parameters class, which carries the model (state) for the computation
"""

import numpy as np


class Parameters:
    """
    The use of this class is to generate a unique input (model) for a dft calculation.
    This includes parameters set by the user (e.g. scf_tol), and any derived quantities
    (e.g. num_planewaves).

    The parameters object is designed to be immutable.
    """

    def __init__(self, **kwargs):

        # Level of approximation used
        self._method = kwargs.get('method', 'dft')

        # Size of real space cell
        self._cell = kwargs.get('cell', 20)
        self._num_planewaves = kwargs.get('num_planewaves', 1000)
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

        # Have v_ext specified by atoms or given explicitly
        self._manual_v_ext = kwargs.get('manual_v_ext', None)

        # Sanity checks
        assert len(self._species) == len(self._positions), 'Each element requires a unique position.'
        assert abs(max(self._positions)) <= self._cell, 'All elements must lie within the primitive unit cell.'
        assert self._num_planewaves % 2 == 0, 'Number of plane-waves must be even.'

        if self._method not in ['h', 'hf', 'dft']:
            raise RuntimeError('Chosen method of {} is not implemented'.format(self._method))

    def smearing_scheme(self, energy):
        """ Ansatz for smearing the occupancies to prevent occupancy-induced instability in SCF iterations """
        return 1 / (np.exp(energy / self._scf_temperature) + 1)

    @property
    def num_planewaves(self):
        return self._num_planewaves

    @property
    def cell(self):
        return self._cell

    @property
    def scf_tol(self):
        return self._scf_tol

    @property
    def scf_history_length(self):
        return self._scf_history_length

    @property
    def scf_step_length(self):
        return self._scf_step_length

    @property
    def num_electrons(self):
        """ Number of electrons (counted s.t. charge neutrality is enforced) """
        num_electrons = 0
        for i in range(len(self._species)):
            num_electrons += self._element_charges[self._species[i]]
        return num_electrons

    @property
    def realspace_grid(self):
        """ Return real space grid for unit cell in Angstroms """
        return np.linspace(-self._cell, self._cell, self._num_planewaves)

    @property
    def planewave_grid(self):
        """ Plane-wave frequences for plane-waves that fit in unit cell: G = 2pi n / R """
        N = int(self._num_planewaves / 2)
        pw_frequencies = [np.pi*n / self._cell
                          for n in
                          range(-N, N)]
        # Ordering defined by numpy's fft
        pw_frequencies[0:N], pw_frequencies[N:2*N] = pw_frequencies[N:2*N], pw_frequencies[0:N]
        return np.asarray(pw_frequencies)

    @property
    def v_ext(self):
        """ External potential: either an atomic potential (Coulomb) or a given functional form """
        if self._manual_v_ext is not None:
            return np.fft.fft(self._manual_v_ext(x=self.realspace_grid))
        else:
            v_ext = np.zeros(self._num_planewaves)
            num_atoms = len(self._species)
            for i in range(num_atoms):
                charge = self._element_charges[self._species[i]]
                v_ext += self.coulomb(charge, self._positions[i])
            return np.fft.fft(v_ext)

    @property
    def k_points(self):
        r""" MP k-point grid within the first BZ: k \in [-pi/a, pi/a] """
        return np.linspace(-np.pi / self._cell, np.pi / self._cell, 1 / self._k_point_spacing)

    def coulomb(self, charge: int, position: float) -> np.ndarray:
        """ Returns the external potential of an ion with a given charge and position """
        assert min(abs(self.realspace_grid - position)) > 0.0001, 'Atom is on a grid-point, please shift'
        return -charge / abs(self.realspace_grid - position)

