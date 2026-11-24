def run(inputs, settings):
    b = inputs["b"]
    psi = inputs["psi"]
    phi_B = inputs["phi_B"]
    f = 0.5 * b * phi_B**2 * (psi - phi_B)**2
    mu_B = b * phi_B * (psi - phi_B) * (psi - 2.0 * phi_B)
    return {"f": f, "mu_B": mu_B}


def symbolic():
    import sympy as sp

    b, psi, phi_B = sp.symbols("b psi phi_B", positive=True)
    f = sp.Rational(1, 2) * b * phi_B**2 * (psi - phi_B)**2
    return {"f": f, "mu_B": sp.diff(f, phi_B)}
