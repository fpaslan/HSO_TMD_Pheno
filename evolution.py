"""
Evolution factor E(bT, Q/Q0), built top to bottom.
HSO approach to Drell-Yan, Aslan et al., arXiv:2401.14266

Version for the E288 fit: input-scale CS kernel, Eq. (74), no RG improvement (Sec. V B).
The RG-improved kernel, Eqs. (77), (80), comes later for the Z0 (Sec. VI).
"""
import numpy as np
import params
from scipy.integrate import quad
from scipy.special import k0

from couplings import alpha_s, CF, Q0

# %% Inputs
###########################    CHANGE BY HAND     ###############################
m_K = 0.3              # [GeV] mass parameter of the CS kernel [Eqs. 70-74]; fixed (Sec. V B)
# fit parameter used here (in params.py): b_K

# ---------------------------------------------------------------- Step 6: evolution factor
# Exponent of Eq. (43):
#   E = exp{ K~(bT; mu_Q0) ln(Q^2/Q0^2)
#            + int_{mu_Q0}^{mu_Q} dmu'/mu' [ 2 gamma(alpha_s(mu'); 1) - ln(Q^2/mu'^2) gamma_K(alpha_s(mu')) ] }
# mu_Q = Q, mu_Q0 = Q0 (Sec. II C, C2 = 1)
def E_evol(bT, Q):
    """Evolution factor E(bT, Q/Q0) from Q0 to Q [Eq. 2; exponent of Eq. 43]."""
    return np.exp(K_tilde(bT) * np.log(Q**2/Q0**2) + mu_integral(Q))


# ---------------------------------------------------------------- 6a: the mu' integral (perturbative, no bT)
#   int_{Q0}^{Q} dmu'/mu' [ 2 gamma - ln(Q^2/mu'^2) gamma_K ]
def mu_integral(Q):
    """Perturbative mu' integral in the exponent of Eq. (43)."""
    def integrand(mu):
        return (2*gamma(mu) - np.log(Q**2/mu**2) * gamma_K(mu)) / mu

    integral, error = quad(integrand, Q0, Q)
    return integral

# Eq. (81): gamma(alpha_s(mu); 1) = 3 CF alpha_s / (2 pi)
def gamma(mu):
    """TMD anomalous dimension gamma(alpha_s(mu); 1) [Eq. 81]."""
    return 3*CF*alpha_s(mu) / (2*np.pi)

# Eq. (68): gamma_K(alpha_s(mu)) = 2 CF alpha_s / pi
def gamma_K(mu):
    """Collins-Soper anomalous dimension gamma_K(alpha_s(mu)) [Eq. 68]."""
    return 2*CF*alpha_s(mu) / np.pi


# ---------------------------------------------------------------- 6b: CS kernel at the input scale
# Eq. (74): K~_inpt(bT; mu_Q0) = 2 pi A_K K0(mK bT) + bK exp(-mK^2 bT^2) + D_K
#   K0 = modified Bessel function of the second kind (text after Eq. 65)
def K_tilde(bT):
    """Collins-Soper kernel at the input scale, K~_inpt(bT; mu_Q0) [Eq. 74]."""
    return 2*np.pi*A_K() * k0(m_K*bT) + params.b_K*np.exp(-m_K**2 * bT**2) + D_K()

# Eq. (71): A_K^(1)(mu_Q0) = alpha_s(mu_Q0) CF / pi^2
def A_K():
    """Coefficient A_K^(1)(mu_Q0) of the K0 term in K~ [Eq. 71]."""
    return alpha_s(Q0)*CF / np.pi**2

# Eq. (73): D_K(mu_Q0) = -bK + (2 alpha_s(mu_Q0) CF / pi) ln(mK / mu_Q0)
def D_K():
    """D_K(mu_Q0): large-bT limit of K~ [Eqs. 73, 76]."""
    return -params.b_K + 2*alpha_s(Q0)*CF/np.pi * np.log(m_K/Q0)






if __name__ == "__main__":
    # check Eq. (76): at large bT, K~ -> D_K
    print("K~(bT=20) =", K_tilde(20.0), "   D_K =", D_K())
    # evolution factor at a few bT, Q = 5.5 GeV
    for bT in [0.1, 1.0, 3.0]:
        print(f"E(bT={bT}, Q=5.5) =", E_evol(bT, 5.5))
