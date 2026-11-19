"""
Input TMD pdf of the proton in bT space, built top to bottom.
HSO approach to Drell-Yan, Aslan et al., arXiv:2401.14266

Version for the E288 fit: input-scale TMD, Eq. (64), no RG improvement (Sec. V B).
The RG-improved TMD, Eq. (78), comes later for the Z0 (Sec. VI).
"""
import numpy as np
import params
from functools import lru_cache
from scipy.integrate import quad
from scipy.special import k0, kv
from couplings import alpha_s, CF, Q0



# %% Inputs
# Values set by hand for this module. Fit parameters (M0, M1, Lam, MX) are in params.py.

core_model = "gauss"     # "gauss" (Eq. 62, Table I) or "spectator" (Eq. 63, Table II)
                         # b_K in params.py goes with the model: 2.12 (Gauss), 1.90 (spectator)

# fixed model parameters
m   = 0.3                # [GeV] mass m_{i,p} = m_{i,p,A} = m_{i,p,B} = m_{i,p,L} [Eq. 50]; Sec. V B
m_g = 0.3                # [GeV] gluon mass m_{g,p} [Eq. 50]; Sec. V B
m_q = 0.0                # [GeV] quark mass in the spectator core [Eq. 63]; set to 0 (Sec. V C)
                      
Mp  = 0.938272           # [GeV] proton mass [Eq. 63]; 

# constants
TF = 0.5                 # QCD color factor


# %% Step 7a: the TMD in bT space
# Eq. (64):

def f_tilde_p(j, x, bT):
    """TMD pdf of the proton at the input scale, f~_inpt,j/p(x, bT; mu_Q0, Q0^2) [Eq. 64].
    j: flavor, x: momentum fraction, bT: transverse position [Eqs. 43, 79]."""
    A, B, Ag, C = coefficients(j, x)
    return (k0(m*bT) * A
            + k0(m*bT) * np.log(Q0**2 * bT / (2*m*np.exp(-gammaE))) * B
            + k0(m_g*bT) * Ag
            + C * f_core(x, bT))

#  inputs
#m   = 0.3                     # [GeV] m_{i,p}: quark mass parameter [Eq. 50]; = m_{i,p,A} = m_{i,p,B} = m_{i,p,L} (Sec. V B)
#m_g = 0.3                     # [GeV] m_{g,p}: gluon mass parameter [Eq. 50] (Sec. V B)
gammaE = np.euler_gamma       # Euler's constant


# %% Step 7b: the core models

def f_core(x, bT):
    """Nonperturbative core in bT space [Eq. 65]; model chosen by core_model."""
    if core_model == "gauss":
        return f_core_gauss(x, bT)
    if core_model == "spectator":
        return f_core_spectator(x, bT)
    raise ValueError(f"unknown core model: {core_model}")


# Gaussian core, Eq. (65): f~_core = exp(-bT^2 MF^2 / 4), with MF -> M0 + M1 log(1/x) (Eq. 89)
def f_core_gauss(x, bT):
    """Gaussian core in bT space [Eqs. 62, 65, 89]."""
    MF = params.M0 + params.M1*np.log(1/x)
    return np.exp(-bT**2 * MF**2 / 4)

         
# Spectator core, Eq. (65):
#   f~_core = (1/4) [(mq + x Mp)^2 - L^2] / [2(mq + x Mp)^2 + L^2] (L bT)^3 K3(L bT)
#           + (3/2) L^2 / [2(mq + x Mp)^2 + L^2] (L bT)^2 K2(L bT)
#   L^2 = (1-x) Lambda^2 + x MX^2 - x(1-x) Mp^2  from (Eq. 63)
#   K2, K3 = modified Bessel functions of the second kind 
def f_core_spectator(x, bT):
    """Spectator-model core in bT space [Eqs. 63, 65]."""
    a  = (m_q + x*Mp)**2
    L2 = (1 - x)*params.Lam**2 + x*params.MX**2 - x*(1 - x)*Mp**2
    z  = np.sqrt(L2) * bT
    return (0.25 * (a - L2) / (2*a + L2) * z**3 * kv(3, z)
            + 1.5 * L2 / (2*a + L2) * z**2 * kv(2, z))

#  inputs: model parameters (Sec. V C)
#m_q = 0.0                     # ["quark" mass in the spectator core (Eq. 63); set to 0 in the paper (Sec. V C)
#Mp  = 0.938272                # [GeV] proton mass 



# %%
@lru_cache(maxsize=None)
def coefficients(j, x):
    """Perturbative-tail coefficients A, B, A^g [Eqs. 51-53] and core coefficient C [Eq. 54]."""
    f = f_col(j, x)

    # Eq. (51): 
    A = alpha_s(Q0)/np.pi * (Pii_conv(j, x) - 3*CF/2 * f)

    # Eq. (52): 
    B = alpha_s(Q0)*CF/np.pi * f

    # Eq. (53): 
    Ag = alpha_s(Q0)/np.pi * Pig_conv(x)

    # Eq. (54), with N_{i/p} = 1 and mu_Q0 = Q0:
    
    C = (f
         - A * np.log(Q0/m)
         - B * np.log(Q0/m) * np.log(Q0**2/(Q0*m))
         - Ag * np.log(Q0/m_g)
         + alpha_s(Q0)/(2*np.pi) * (CDii_conv(j, x) + CDig_conv(x)))
    return A, B, Ag, C


# %% Step 7c ingredients: convolutions

def Pii_conv(j, x):
    """Convolution (P_ii x f_j)(x) [Eqs. 51, 55, 61]"""
    f = f_col(j, x)

    def integrand(z):
        return ((1 + z**2) * f_col(j, x/z) / z - 2*f) / (1 - z)

    integral, error = quad(integrand, x, 1)
    return CF * (integral + 2*f*np.log(1 - x) + 3/2*f)


def Pig_conv(x):
    """Convolution (P_ig x f_g)(x) [Eqs. 53, 56, 61]"""
    def integrand(z):
        return TF * (z**2 + (1 - z)**2) * f_col("g", x/z) / z

    integral, error = quad(integrand, x, 1)
    return integral


def CDii_conv(j, x):
    """Convolution (C_Delta^{i/i} x f_j)(x) [Eqs. 54, 57, 61]"""
    def integrand(z):
        return CF * (1 - z) * f_col(j, x/z) / z

    integral, error = quad(integrand, x, 1)
    return integral - CF * np.pi**2/12 * f_col(j, x)


def CDig_conv(x):
    """Convolution (C_Delta^{i/g} x f_g)(x) [Eqs. 54, 58, 61]"""
    def integrand(z):
        return 2*TF  * (1 - z) * f_col("g", x/z) 

    integral, error = quad(integrand, x, 1)
    return integral


# %% bottom layer: collinear pdfs
# f_{i/p}(x; mu_Q0): MMHT2014 NLO (Sec. V), at mu_Q0 = Q0.  j = "u", "ubar", ..., "g"
#For the MS collinear pdfs, we use the NLO extraction MMHT2014 of Ref. [45], accessible through LHAPDF6 [43].

def f_col(j, x):
    """Collinear pdf f_{j/p}(x; mu_Q0), MMHT2014 NLO [Eq. 4; Sec. V]."""
    if j == "g":
        return (1 - x)**5     # placeholder shape (not physical); will come from pdfs.py
    else:
        return (1 - x)**3         # placeholder shape (not physical); will come from pdfs.py





# %%

if __name__ == "__main__":
    j, x = "u", 0.3
    A, B, Ag, C = coefficients(j, x)
    print(f"A = {A:.4f}  B = {B:.4f}  A^g = {Ag:.4f}  C = {C:.4f}")
    for bT in [0.05, 0.5, 2.0, 5.0]:
        print(f"f~(bT={bT}) =", f_tilde_p(j, x, bT))
        
        
        
        
        
        
        
        
        