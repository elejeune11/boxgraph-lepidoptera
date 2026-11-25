import math

def run(inputs, settings):
    tau = inputs["eta_rot"] / inputs["K_S"]
    frac = 1.0 - math.exp(-inputs["t"] / tau)
    return {
        "S_t": inputs["S_0"] * frac,
        "S_over_S0": frac,
        "tau_relax": tau,
    }

def symbolic():
    import sympy as sp
    S_0, K_S, eta_rot, t = sp.symbols("S_0 K_S eta_rot t", positive=True)
    tau = eta_rot / K_S
    frac = 1 - sp.exp(-t / tau)
    return {"S_t": S_0 * frac, "S_over_S0": frac, "tau_relax": tau}
