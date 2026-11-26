import math


def run(inputs, settings):
    f_comp = float(inputs["f_comp"])  # pN
    a = float(inputs["a"])            # nm
    L = float(inputs["L"])            # um
    # Convert L from um to nm (L * 1e3) so a^2 / L_nm has units of nm
    Delta_E_comp = -(f_comp / 16.0) * (a ** 2) * (math.pi ** 2 / (L * 1e3))
    return {"Delta_E_comp": Delta_E_comp}


def symbolic():
    import sympy as sp

    f_comp, a, L = sp.symbols("f_comp a L", positive=True)
    return {"Delta_E_comp": -(f_comp / 16) * a**2 * (sp.pi**2 / (L * 1000))}
