import numpy as np

"""
Input to the calculation in the form of a parameters class
"""


class parameters(object):
    def __init__(self,*args,**kwargs):

        # Level of approximation used
        self.method = 'dft'

        # Size of real space cell
        self.cell = 20
        self.num_planewaves = 101
        self.kpoint_spacing = 0.05

        # Real space grid
        self.grid = np.linspace(-0.5*self.cell, 0.5*self.cell, self.num_planewaves)

        def get_planewave_frequencies():
            print('hello')

        def get_kpoints():
            print('hello')

        # List of species + position
        self.species = ['Li']
        self.position = [0]

        # Define each element with corresponding atomic number
        self.element_charges = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5,
                                'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10}

        # SCF
        self.scf_tol = 1e-10
        self.history_length = 10
        self.step_length = 1

        # Have v_ext specified by atoms or given explicitly
        self.manual_v_ext = False
        self.v_ext = 0.5*(0.25**2)*self.grid**2

        # Number of atoms in the calculation
        self.num_atoms = len(self.species)

        # Number of `electrons'
        self.num_particles = 0
        for i in range(0,len(self.species)):
            self.num_particles += self.element_charges[self.species[i]]



