def run(inputs, settings):
    k = inputs["k"]
    rate_0 = inputs["gamma_I"] / (inputs["eta"] * inputs["R_0"])
    return {"Omega": rate_0 * (k**2) * (1.0 - (1.0 + 2.0 * inputs["omega"]) * (k**2))}


def symbolic():
    import sympy as sp

    gamma_I, eta, R_0, omega, k = sp.symbols("gamma_I eta R_0 omega k", positive=True)
    return {"Omega": (gamma_I / (eta * R_0)) * k**2 * (1 - (1 + 2 * omega) * k**2)}
