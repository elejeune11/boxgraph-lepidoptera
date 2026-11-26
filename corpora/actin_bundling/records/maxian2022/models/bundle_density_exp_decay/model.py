import math


def run(inputs, settings):
    return {"rho_B": inputs["A"] * math.exp(-(inputs["t"] - inputs["t0"]) / inputs["tau_dec"])}


def symbolic():
    import sympy as sp

    t, A, t0, tau_dec = sp.symbols("t A t0 tau_dec", positive=True)
    return {"rho_B": A * sp.exp(-(t - t0) / tau_dec)}
