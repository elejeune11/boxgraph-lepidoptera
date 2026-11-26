"""Intermediate shear-dominated regime scaling law (p. 5 col 1): kappa_B propto N * k_x * L^2.

Printed without a prefactor (prefactor: unknown); exact Eq. (2) asymptotic limit is
kappa_B = (1 / (12*c)) * N * (d_f + t)^2 * (k_x / delta) * L^2.
"""


def run(inputs, settings):
    return {"kappa_B": inputs["N"] * inputs["k_x"] * (inputs["L"] ** 2)}


def symbolic():
    import sympy as sp

    N, k_x, L = sp.symbols("N k_x L", positive=True)
    return {"kappa_B": N * k_x * L**2}
