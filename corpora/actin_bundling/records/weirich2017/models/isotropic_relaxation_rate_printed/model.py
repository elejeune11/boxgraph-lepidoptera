def run(inputs, settings):
    # As printed on SI p. 8 col 2 (dimensionally inverted: r*eta/gamma_I has units of s,
    # and gamma_I/eta has units of m/s rather than s). Evaluated verbatim in SI for (r, eta, gamma_I).
    r_m = inputs["r"] * 1e-6
    gamma_Nm = inputs["gamma_I"] * 1e-9
    dL_dt = -inputs["L"] * (r_m * inputs["eta"] / gamma_Nm)
    tau = gamma_Nm / inputs["eta"]
    return {"dL_dt": dL_dt, "tau": tau}


def symbolic():
    import sympy as sp

    L, r, eta, gamma_I = sp.symbols("L r eta gamma_I", positive=True)
    # r_m / gamma_Nm = (r * 1e-6) / (gamma_I * 1e-9) = 1000 * r / gamma_I
    return {
        "dL_dt": -sp.Integer(1000) * L * r * eta / gamma_I,
        "tau": sp.Rational(1, 10**9) * gamma_I / eta,
    }
