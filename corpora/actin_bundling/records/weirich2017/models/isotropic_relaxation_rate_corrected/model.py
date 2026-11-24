def run(inputs, settings):
    # gamma_I [nN/m = 1e-9 N/m] / (r [um = 1e-6 m] * eta [Pa*s]) = 1e-3 * gamma_I / (r * eta) [1/s]
    rate_per_s = 1e-3 * inputs["gamma_I"] / (inputs["r"] * inputs["eta"])
    return {
        "dL_dt": -inputs["L"] * rate_per_s,
        "tau": 1.0 / rate_per_s,
    }


def symbolic():
    import sympy as sp

    L, r, eta, gamma_I = sp.symbols("L r eta gamma_I", positive=True)
    return {
        "dL_dt": -sp.Rational(1, 1000) * L * gamma_I / (r * eta),
        "tau": sp.Integer(1000) * r * eta / gamma_I,
    }
