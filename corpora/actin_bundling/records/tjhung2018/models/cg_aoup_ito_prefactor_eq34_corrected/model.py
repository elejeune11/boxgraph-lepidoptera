def run(inputs, settings):
    # Corrected Stratonovich-to-Ito drift: s = +1 (+tau*v*grad_i v) -> prefactor 1 - 1/2 = 1/2 (Eq. 39).
    s = 1.0
    return {"log_prefactor": 1.0 - 0.5 * s}


def symbolic():
    import sympy as sp

    tau, v0 = sp.symbols("tau v0", positive=True)
    return {"log_prefactor": sp.Rational(1, 2) + 0 * tau * v0}
