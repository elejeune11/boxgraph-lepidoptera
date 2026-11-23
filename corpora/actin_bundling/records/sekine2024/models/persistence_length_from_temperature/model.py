"""Composite bending persistence length substituting k_B T = k_B * T (Table S1)."""

K_B_PNUM_PER_K = 1.380649e-5


def run(inputs, settings):
    kBT = K_B_PNUM_PER_K * inputs["T"]
    return {"l_bend": inputs["k_bend"] * inputs["L0"] / (2.0 * kBT)}


def symbolic():
    import sympy as sp

    k_bend, L0, T = sp.symbols("k_bend L0 T", positive=True)
    return {"l_bend": k_bend * L0 / (2 * sp.Float(K_B_PNUM_PER_K) * T)}
