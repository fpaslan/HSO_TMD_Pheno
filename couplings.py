"""
Couplings and constants shared by all modules.
HSO approach to Drell-Yan, Aslan et al., arXiv:2401.14266
"""


# %% Couplings
# electromagnetic coupling
def alpha_em(Q):
    """Electromagnetic coupling """
    return 1/137.036          # placeholder


# strong coupling
def alpha_s(mu):
    """Strong coupling"""
    return 0.2                # placeholder


# %% Constants

Nc = 3                        # number of colors
CF = 4/3                      # QCD color factor 
Q0 = 4.0                      # [GeV] input scale of the TMDs (Sec. V B)


# %% Charges
charges = {"u": 2/3, "d": -1/3, "s": -1/3, "c": 2/3, "b": -1/3}   