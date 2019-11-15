""" Example code for 1D Silicon """

# Create Parameters object (see Documentation) with various inputs
# params = Parameters(species=['Si'],
#                      positions=[0],
#                      num_planewaves=1000)

# Default constuction of density is overlapping of atomic oribtals
# input ndarray to construct
# density = Density() 

# Create abstract Hamiltonian
# We create the gamma-point hamiltonian for now.
# hamiltonian = Hamiltonian(density=density
#                           k_point=0)

# Compute the gamma point wavefunction
# wavefunction = hamiltonian.get_eigenfunctions()

# And its density
# density = wavefunction.get_density()

# non self-consistent energy of 1D silicon
# energy = hamiltonian.total_energy(wavefunctions)

# iterate toward SCF
# for it in scf

# Out is our scf density, that we can input to hamiltonian to get various thigns: wvfns, energies, etc.

