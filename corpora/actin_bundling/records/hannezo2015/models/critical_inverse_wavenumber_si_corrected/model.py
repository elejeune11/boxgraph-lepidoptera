import math


def run(inputs, settings):
    D_si = inputs["D"] * 1e-12
    ratio = D_si * inputs["tau"] * inputs["xi"] / inputs["eta"]
    inv_qc = inputs["L"] * (ratio ** 0.25)
    return {"inv_qc": inv_qc, "q_c": 1.0 / inv_qc}


def symbolic():
    import sympy as sp

    L, D, tau, xi, eta = sp.symbols("L D tau xi eta", positive=True)
    D_si = D * sp.Rational(1, 10**12)
    inv_qc = L * (D_si * tau * xi / eta) ** sp.Rational(1, 4)
    return {"inv_qc": inv_qc, "q_c": 1 / inv_qc}
