"""
E288 cross section, Eq. (86), built top to bottom.
HSO approach to Drell-Yan, Aslan et al., arXiv:2401.14266
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import j0

from couplings import alpha_em, alpha_s, Nc, CF, Q0, charges
from switches import placeholder                    # placeholder switches, all modules
from evolution import E_evol as E_evol_module      # real code, Step 6
from hadrons import f_tilde                         # TMDs of the hadrons: proton, copper target (hadrons.py)






# %% Inputs
# Experimental setup: 
    
###########################    CHANGE BY HAND     ###############################

Ebeam = 400.0                 # [GeV] beam energy: 200, 300 or 400 for E288 (Sec. V A)
h_a = "p"                     # hadron a: beam proton
h_b = "Cu"                    # hadron b: copper target nucleon


bT_max = 30.0                 # [GeV^-1] upper limit of the bT integral


#Flavors
flavors = ["u", "d", "s"]     # flavors summed in Eq. (43)



# %% Kinematics for E288
# Quantities fixed by the experiment, one value per beam energy.
# E288 values, from  arXiv:1901.06916: sqrt(s) in the E288 section, <y_h> in Figs. 6-9 (called eta there)

sqrt_s_table = {200.0: 19.4, 300.0: 23.8, 400.0: 27.4}     # [GeV]  c.m. energy
yh_table     = {200.0: 0.40, 300.0: 0.21, 400.0: 0.03}     # <y_h> (Eq. 85)

s  = sqrt_s_table[Ebeam]**2     
yh = yh_table[Ebeam]            

# momentum fractions, from Eq. (15)
def x_ab(Q, qT):
    """Momentum fractions x_a, x_b [Eqs. 10, 15]."""
    if placeholder["x_ab"] is not None:
        return placeholder["x_ab"]
    root = np.sqrt(Q**2 / (s * (1 + qT**2/Q**2)))
    return root*np.exp(yh), root*np.exp(-yh)

# %% 





























# %% Step 1: observable
# Eq. (86): (1/pi) * integral over the Q bin, at fixed <qT> and <y_h>
def O_E288(Qmin, Qmax, qT):
    """E288 observable: (1/pi) int dQ^2 (dsigma) over a Q bin [Eq. 86].
    Qmin,Qmax, <qT> and <y_h> is given by the experimental collaboration [Eq.85] """
    integral, error = quad(dsigma, Qmin**2, Qmax**2, args=(qT,))
    return integral / np.pi


# %% Step 2: integrand
# Eq. (29), with F2 = 0 (Eq. 37)
def dsigma(Q2, qT):
    """Cross section in [Eq. 29] with Q2 = Q^2"""
    if placeholder["dsigma"] is not None:
        return placeholder["dsigma"]
    Q = np.sqrt(Q2)
    return 2*np.pi**2 * alpha_em(Q)**2 / (3*s*Q2) * 2*F1(Q, qT)

#  inputs
#  alpha_em(Q): from couplings.py Has a placeholder value for now. Fix later!!!
#  s: c.m. energy squared -- see Kinematics cell at the very top


# %% Step 3: structure function
# Eq. (43): sum over flavors j 
#   F1 = sum_j e_j^2 |H|^2 / (4 pi^2 Nc) * [bT integral for flavor j]
def F1(Q, qT):
    """Structure function F^1_UU  [Eq 43]."""
    if placeholder["F1"] is not None:
        return placeholder["F1"]
    total = 0.0
    for j in flavors:
        total += charges[j]**2 * bT_integral(j, Q, qT)
    return H2(Q) / (4*np.pi**2*Nc) * total



#  inputs
#  Nc: from couplings.py
#  charges e_j
#  which flavors are summed  -> placeholder list, we can modify 



# %% Step 4: hard factor
# Eq. (34) truncated at O(alpha_s), with H^(1) from Eq. (35); mu_Q = Q (Sec. II C, C2 = 1)
def H2(Q):
    """Hard factor |H_jjbar|^2 [Eqs. 34-35]."""
    if placeholder["H2"] is not None:
        return placeholder["H2"]
    return 1 + alpha_s(Q)/(4*np.pi) * CF * (-16 + 7*np.pi**2/3)

#  inputs
#  alpha_s(mu), CF: from couplings.py


# %% Step 5: bT integral
# The integral in Eq. (43) for one flavor j:
#   int d^2bT e^{i qT.bT} [ f~_{j/ha}(xa) f~_{jbar/hb}(xb) + (a <-> b) ] * E(bT, Q/Q0)
#
#   int d^2bT e^{i qT.bT} (...) = 2 pi int_0^inf dbT bT J0(qT bT) (...)
def bT_integral(j, Q, qT):
    """bT integral in Eq. (43)."""
    if placeholder["bT_integral"] is not None:
        return placeholder["bT_integral"]
    xa, xb = x_ab(Q, qT)
    jbar = antiquark[j]          # antiquark of flavor j, e.g. "u" -> "ubar"

    def integrand(bT):
        tmds = (f_tilde(j, h_a, xa, bT) * f_tilde(jbar, h_b, xb, bT)     # f~_{j/ha} f~_{jbar/hb}
              + f_tilde(j, h_b, xb, bT) * f_tilde(jbar, h_a, xa, bT))    # (a <-> b)
        return bT * j0(qT*bT) * tmds * E_evol(bT, Q)

 
    integral, error = quad(integrand, 0.0, bT_max, limit=200)
    #print(f"bT integral, flavor {j}: {integral:.4e} +- {error:.1e}")
    return 2*np.pi * integral
# limit = max number of subintervals quad may use (default 50).
# J0(qT bT) oscillates, so quad needs many short subintervals to follow the wiggles;
# with only 50 it may stop before reaching its accuracy. 200 lets it finish.
    integral, error = quad(integrand, 0.0, bT_max, limit=200)
#  inputs
#bT_max = 30.0                 # [GeV^-1] upper limit of the bT integral

# antiquark label for each quark flavor
antiquark = {"u": "ubar", "d": "dbar", "s": "sbar"}




# %% Step 5 ingredient: evolution factor
# E(bT, Q/Q0): Eq. (2), the exponential in Eq. (43) -- real code in evolution.py (Step 6)
def E_evol(bT, Q):
    """Evolution factor E(bT, Q/Q0) from Q0 to Q [Eq. 2; exponent of Eq. 43]."""
    if placeholder["E_evol"] is not None:
        return placeholder["E_evol"]
    return E_evol_module(bT, Q)



# %% Run # O_E288(Qmin, Qmax, qT):
if __name__ == "__main__":
    # Qmin = 5, Qmax = 6, qT = 1.0, y = 0.03: placeholders, not E288 values 
    print(O_E288(5, 6, qT=1.0))
