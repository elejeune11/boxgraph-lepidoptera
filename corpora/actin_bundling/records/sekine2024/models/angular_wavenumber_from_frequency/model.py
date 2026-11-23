"""Convert spatial frequency mu/(2*pi) (Supp. Table 4) to angular wavenumber mu (Methods p. 11)."""

import math


def run(inputs, settings):
    return {"mu": 2.0 * math.pi * inputs["mu_over_2pi"]}


def symbolic():
    import sympy as sp

    mu_over_2pi = sp.Symbol("mu_over_2pi", positive=True)
    return {"mu": 2 * sp.pi * mu_over_2pi}
