import math


def run(inputs, settings):
    kBT = float(inputs["kBT"])  # pN*nm
    L_p = float(inputs["L_p"])  # um
    a = float(inputs["a"])      # nm
    L = float(inputs["L"])      # um
    # Convert a from nm to um (a * 1e-3) so L_p * a_um^2 / L^3 is dimensionless
    a_um = a * 1e-3
    E_bend = (kBT * L_p / 64.0) * (a_um ** 2) * (math.pi ** 4 / L ** 3)
    return {"E_bend": E_bend}


def symbolic():
    import sympy as sp

    kBT, L_p, a, L = sp.symbols("kBT L_p a L", positive=True)
    a_um = a * sp.Rational(1, 1000)
    return {"E_bend": (kBT * L_p / 64) * a_um**2 * (sp.pi**4 / L**3)}
