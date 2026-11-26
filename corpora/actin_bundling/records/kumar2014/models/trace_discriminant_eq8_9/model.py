def run(inputs, settings):
    D = inputs["D"]
    k = inputs["k"]
    ell = inputs["ell"]
    alpha = inputs["alpha"]
    P = inputs["P"]
    nu = inputs["nu"]
    a = inputs["A0_fA"]
    b = inputs["I0_fI"]
    Pi = P / (1.0 + nu * (k * ell) ** 2)
    tr_L = -D * (k**2) * ((1.0 + alpha) - Pi * (a + b))
    disc_L = (D**2) * (k**4) * (
        (1.0 - alpha) ** 2
        + (Pi**2) * (a + b) ** 2
        - 2.0 * Pi * (1.0 - alpha) * (a - b)
    )
    tau = (ell**2) / D
    return {
        "tr_L": tr_L,
        "disc_L": disc_L,
        "tr_tau_L": tau * tr_L,
        "disc_tau_L": (tau**2) * disc_L,
    }


def symbolic():
    import sympy as sp

    D, k, ell, alpha, P, nu = sp.symbols("D k ell alpha P nu", positive=True)
    A0_fA, I0_fI = sp.symbols("A0_fA I0_fI", real=True)
    Pi = P / (1 + nu * k**2 * ell**2)
    tr_L = -D * k**2 * ((1 + alpha) - Pi * (A0_fA + I0_fI))
    disc_L = D**2 * k**4 * (
        (1 - alpha) ** 2
        + Pi**2 * (A0_fA + I0_fI) ** 2
        - 2 * Pi * (1 - alpha) * (A0_fA - I0_fI)
    )
    tau = ell**2 / D
    return {
        "tr_L": tr_L,
        "disc_L": disc_L,
        "tr_tau_L": sp.simplify(tau * tr_L),
        "disc_tau_L": sp.simplify(tau**2 * disc_L),
    }
