import math


def run(inputs, settings):
    kBT = float(inputs["kBT"])  # pN*nm
    L_p = float(inputs["L_p"])  # um
    a = float(inputs["a"])      # nm
    L = float(inputs["L"])      # um

    f_buckle_pN = (math.pi ** 2) * kBT * L_p / (4.0 * L ** 2) * 1e-3
    Delta_E_tot = (math.pi ** 2 * a ** 2 / (16.0 * L * 1e3)) * f_buckle_pN
    return {"Delta_E_tot": Delta_E_tot}


def symbolic():
    import sympy as sp

    kBT, L_p, a, L = sp.symbols("kBT L_p a L", positive=True)
    f_buckle_pN = sp.pi**2 * kBT * L_p / (4 * L**2) * sp.Rational(1, 1000)
    return {"Delta_E_tot": (sp.pi**2 * a**2 / (16 * L * 1000)) * f_buckle_pN}
