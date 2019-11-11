"""
Class designed to keep and update the relevant state (e.g. density history)
of a calculation en route toward self-consistency
"""

import numpy as np


class SCF:
    """
    self stores state, namely, history of densities and residuals for quasi-Newton methods
    """
    def __init__(self, params):
        self._history_length = params.scf_history_length

    def __next__(self, params):
        return None
        #self._density_out = kohn_sham_map(params, self._density)
        # Add to list of in/out densities
        # return self.density = pulay()

    @staticmethod
    def kohn_sham_map(params, density_in):
        # return density_out
        return None

    def pulay_update(self):
        # Compute rho_out
        return None