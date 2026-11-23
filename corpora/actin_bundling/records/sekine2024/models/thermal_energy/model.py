"""Thermal energy scale k_B T in pN*um (SI Notes Table S1)."""

# CODATA 2018 exact Boltzmann constant: 1.380649e-23 J/K = 1.380649e-5 pN*um/K
K_B_PNUM_PER_K = 1.380649e-5


def run(inputs, settings):
    return {"kBT": K_B_PNUM_PER_K * inputs["T"]}


def symbolic():
    import sympy as sp

    T = sp.Symbol("T", positive=True)
    return {"kBT": sp.Float(K_B_PNUM_PER_K) * T}
