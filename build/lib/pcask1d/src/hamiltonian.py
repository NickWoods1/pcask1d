import numpy as np
import matplotlib.pyplot as plt

class hamiltonian(object):

    # Constructor function
    def __init__(self, params, *args, **kwargs):

        # Default gamma point
        self.kpoint = kwargs.pop('kpoint',0)
        self.representation = np.zeros((params.num_planewaves,params.num_planewaves))

        planewave_frequencies = params.get_planewave_frequencies()

        # Kinetic energy
        for i in range(params.num_planewaves-1):
            self.representation[i,i] += 1 / abs(planewave_frequencies[i] + self.kpoint)**2

        # External potential
        x = np.linspace(-10,10,100)
        """
        v_ext = np.fft.fft(np.exp(-x**2))
        freq = np.fft.fftfreq(x.shape[-1])
        print(freq)
        plt.plot(np.exp(-x**2))
        #plt.plot(v_ext*np.exp(1.0j*freq*x))

        a = np.fft.ifft(v_ext)
        plt.plot(a)
        plt.show()

        plt.plot(freq, v_ext.real)
        plt.plot(freq, v_ext.imag)
        plt.show()
        """

    def set_hartree(self, density):
        print('hello')

    def set_vxc(self, density):
        print('hello')

    def set_kpoint(self, kpoint):
        self.kpoint = kpoint

    def get_eigenfunctions(self):
        print('hello')