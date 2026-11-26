import math


def run(inputs, settings):
    # Verbatim printed Fig. S5 legend: 1.4 * exp((-t - 30)/60) and 1.4 * exp((-t - 10)/25)
    return {"rho_B": inputs["A"] * math.exp((-inputs["t"] - inputs["t0"]) / inputs["tau_dec"])}


def symbolic():
    import sympy as sp

    t, A, t0, tau_dec = sp.symbols("t A t0 tau_dec", positive=True)
    return {"rho_B": A * sp.exp((-t - t0) / tau_dec)}
