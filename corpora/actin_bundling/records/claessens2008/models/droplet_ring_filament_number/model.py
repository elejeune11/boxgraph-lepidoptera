import math

N_A = 6.02214076e23  # 1/mol (exact SI)


def run(inputs, settings):
    # c_a [mg/ml] = 1e-15 g/um^3; V = (pi/6)*D_d^3 [um^3]; d [nm] = 1e-3 um
    # Total factor: 1e-15 * 1e-3 = 1e-18
    c_a = inputs["c_a"]
    D_d = inputs["D_d"]
    M_actin = inputs["M_actin"]
    d = inputs["d"]
    L = (math.pi / 6.0) * (D_d**3) * (c_a / M_actin) * N_A * d * 1e-18
    n_f = L / (math.pi * D_d)
    return {"L": L, "n_f": n_f}


def symbolic():
    import sympy as sp

    c_a, D_d, M_actin, d = sp.symbols("c_a D_d M_actin d", positive=True)
    coeff = sp.Float(N_A * 1e-18)
    L = (sp.pi / 6) * D_d**3 * (c_a / M_actin) * coeff * d
    n_f = sp.Rational(1, 6) * D_d**2 * (c_a / M_actin) * coeff * d
    return {"L": L, "n_f": n_f}
