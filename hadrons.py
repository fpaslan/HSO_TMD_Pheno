"""
TMD pdfs of the hadrons in the collision: proton, neutron, copper target.
HSO approach to Drell-Yan, Aslan et al., arXiv:2401.14266
"""
from switches import placeholder                    # placeholder switches, all modules
from tmd import f_tilde_p as f_tilde_p_module       # real code, Step 7


# %% Hadron TMDs
# Input TMD f~_{j/h}(x, bT; mu_Q0, Q0^2) for hadron h = "p" or "Cu"
# Z, A for each nuclear target [Eq. 88]
nuclei = {"Cu": (29, 63)}        # copper (Sec. V A); add others here by hand 

def f_tilde(j, h, x, bT):
    if placeholder["f_tilde"] is not None:
        return placeholder["f_tilde"]
    if h == "p":
        return f_tilde_p(j, x, bT)
    if h in nuclei:
        Z, A = nuclei[h]
        return f_tilde_target(j, x, bT, Z, A)
    raise ValueError(f"unknown hadron: {h}")
    
    
    
# %%
def f_tilde_target(j, x, bT, Z, A):
    """TMD pdf of a target[Eq. 88]."""
    if placeholder["f_tilde_target"] is not None:
        return placeholder["f_tilde_target"]
    return Z/A * f_tilde_p(j, x, bT) + (A - Z)/A * f_tilde_n(j, x, bT)
 

# %%
# Standard isospin: u <-> d, ubar <-> dbar, s unchanged.
def f_tilde_n(j, x, bT):
    """TMD pdf of the neutron, from the proton by isospin [Eq. 88; refs. 46, 47]."""
    if placeholder["f_tilde_n"] is not None:
        return placeholder["f_tilde_n"]
    swap = {"u": "d", "d": "u", "ubar": "dbar", "dbar": "ubar"}
    return f_tilde_p(swap.get(j, j), x, bT)
# swap: isospin table, proton flavor <-> neutron flavor (u <-> d, ubar <-> dbar)
# swap.get(j, j): look up j in swap; if j is not in the table (s, sbar, g), keep j unchanged, second j is the fallback

# %%
# Proton TMD in bT space: Eq. (64) -- call from tmd.py (Step 7)
def f_tilde_p(j, x, bT):
    """TMD pdf of the proton at the input scale [Eq. 64]."""
    if placeholder["f_tilde_p"] is not None:
        return placeholder["f_tilde_p"]
    return f_tilde_p_module(j, x, bT)
