"""Connector bending persistence length l^{X,bend} = k^{X,bend} * L_0^X / (2 * k_B T) (Table S1)."""


def run(inputs, settings):
    return {"l_bend": inputs["k_bend"] * inputs["L0"] / (2.0 * inputs["kBT"])}


def symbolic():
    import sympy as sp

    k_bend, L0, kBT = sp.symbols("k_bend L0 kBT", positive=True)
    return {"l_bend": k_bend * L0 / (2 * kBT)}
