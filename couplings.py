"""
Couplings and constants shared by all modules.
HSO approach to Drell-Yan, Aslan et al., arXiv:2401.14266
"""


# %% Couplings
# electromagnetic coupling, used in Eq. (29)
def alpha_em(Q):
    """Electromagnetic coupling [used in Eq. 29]; no equation in paper."""
    return 1/137.036          # placeholder; no equation in paper (alpha_em(MZ) = 1/127.951 quoted after Eq. 39)


# strong coupling, used in Eqs. (34)-(35) and in the evolution
def alpha_s(mu):
    """Strong coupling [Sec. V, MMHT2014]; no equation in paper."""
    return 0.2                # placeholder; from MMHT2014 NLO set (Sec. V)


# %% Constants

Nc = 3                        # number of colors [Eqs. 36, 43]; value not in paper
CF = 4/3                      # QCD color factor [Eqs. 35, 55, 57, 68, 81]; value not in paper
Q0 = 4.0                      # [GeV] input scale of the TMDs (Sec. V B)


# %% Charges
charges = {"u": 2/3, "d": -1/3, "s": -1/3, "c": 2/3, "b": -1/3}    # quark charges e_j; values not in paper