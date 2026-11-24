def run(inputs, settings):
    a = inputs["a"]
    b = inputs["b"]
    psi = inputs["psi"]
    phi_A = inputs["phi_A"]
    phi_B = inputs["phi_B"]
    f = 0.5 * a * phi_A**2 + (b / (2.0 * psi**2)) * phi_B**2 * (psi - phi_B)**2
    mu_A = a * phi_A
    mu_B = (b / psi**2) * phi_B * (psi - phi_B) * (psi - 2.0 * phi_B)
    return {"f": f, "mu_A": mu_A, "mu_B": mu_B}


def symbolic():
    import sympy as sp

    a, b, psi, phi_A, phi_B = sp.symbols("a b psi phi_A phi_B", positive=True)
    f = sp.Rational(1, 2) * a * phi_A**2 + b / (2 * psi**2) * phi_B**2 * (psi - phi_B)**2
    return {
        "f": f,
        "mu_A": sp.diff(f, phi_A),
        "mu_B": sp.diff(f, phi_B),
    }
