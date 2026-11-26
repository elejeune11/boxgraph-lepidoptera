import math


def run(inputs, settings):
    kBT = float(inputs["kBT"])        # pN*nm
    L_p = float(inputs["L_p"])        # um
    a = float(inputs["a"])            # nm
    L = float(inputs["L"])            # um
    f_comp = float(inputs["f_comp"])  # pN

    # f_buckle in pN: (pi^2 * kBT * L_p / (4 * L^2)) in fN -> * 1e-3 for pN
    f_buckle_pN = (math.pi ** 2) * kBT * L_p / (4.0 * L ** 2) * 1e-3
    # a^2 / L_nm = a^2 / (L * 1e3) in nm
    Delta_E_tot = (a ** 2) * (math.pi ** 2 / (16.0 * L * 1e3)) * (f_buckle_pN - f_comp)
    return {"Delta_E_tot": Delta_E_tot}


def symbolic():
    import sympy as sp

    kBT, L_p, a, L, f_comp = sp.symbols("kBT L_p a L f_comp", positive=True)
    f_buckle_pN = sp.pi**2 * kBT * L_p / (4 * L**2) * sp.Rational(1, 1000)
    return {"Delta_E_tot": a**2 * (sp.pi**2 / (16 * L * 1000)) * (f_buckle_pN - f_comp)}
