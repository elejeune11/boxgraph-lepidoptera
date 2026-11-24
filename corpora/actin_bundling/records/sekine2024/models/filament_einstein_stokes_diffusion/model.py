"""Single-filament Einstein-Stokes center-of-mass diffusivity from Sekine's cylindrical drag (SI Notes I & Table S1).

Combines thermal_energy (k_B T = K_B * T in pN*um) with two-bead total cylindrical friction
2 * gamma^f = 6 * pi * eta_cyto * (3 * d_f + 2 * L0_f / 5) [pN*s/um] to yield
D_fil = k_B T / (2 * gamma^f) [um^2/s].
"""

import math

# CODATA 2018 exact Boltzmann constant: 1.380649e-23 J/K = 1.380649e-5 pN*um/K
K_B_PNUM_PER_K = 1.380649e-5


def run(inputs, settings):
    kBT = K_B_PNUM_PER_K * inputs["T"]
    gamma_tot = 6.0 * math.pi * inputs["eta_cyto"] * (3.0 * inputs["d_f"] + 0.4 * inputs["L0_f"])
    return {"D_fil": kBT / gamma_tot}


def symbolic():
    import sympy as sp

    T, eta_cyto, d_f, L0_f = sp.symbols("T eta_cyto d_f L0_f", positive=True)
    gamma_tot = 6 * sp.pi * eta_cyto * (3 * d_f + sp.Rational(2, 5) * L0_f)
    return {"D_fil": sp.Float(K_B_PNUM_PER_K) * T / gamma_tot}
