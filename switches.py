"""
Placeholder switches shared by all modules.
HSO approach to Drell-Yan, Aslan et al., arXiv:2401.14266

number (or (xa, xb) for x_ab) = closed: the function returns that value, nothing below it runs
None                           = open:   the function uses the real code
"""

placeholder = {
    # e288_xsec.py
    "dsigma":         None,
    "F1":             None,
    "H2":             None,
    "bT_integral":    None,
    "x_ab":           None,     # e.g. (0.2, 0.2)
    "E_evol":         None,
    # hadrons.py
    "f_tilde":        None,
    "f_tilde_target": None,
    "f_tilde_n":      None,
    "f_tilde_p":      None,
}
