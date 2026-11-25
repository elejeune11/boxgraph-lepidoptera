import math


def _F_single(x, W):
    return (x ** 2) - W * (1.0 - 1.0 / math.sqrt(1.0 + x ** 2))


def run(inputs, settings):
    W = inputs["W"]
    V_ratio = inputs["V_ratio"]
    f_V = inputs["f_V"]
    gamma_pN_um = inputs["gamma"] * 1e-3  # nN/m -> pN/um
    a = inputs["a"]                       # um

    x1 = (max(0.0, f_V) * V_ratio) ** (1.0 / 3.0)
    x2 = (max(0.0, 1.0 - f_V) * V_ratio) ** (1.0 / 3.0)
    F_tot_dimless = _F_single(x1, W) + _F_single(x2, W)
    F_tot = gamma_pN_um * (a ** 2) * F_tot_dimless
    u_L = x1 / (x1 + x2) if (x1 + x2) > 0.0 else 0.0

    # Symmetric second derivative d^2(F_tot/(gamma*a^2))/d(R_1/a)^2 at x1 = x2 = x*
    x_s = (0.5 * V_ratio) ** (1.0 / 3.0)
    d2F_sym = -4.0 + 2.0 * W * (1.0 + 4.0 * x_s ** 2) / ((1.0 + x_s ** 2) ** 2.5)
    return {
        "F_tot_dimless": F_tot_dimless,
        "F_tot": F_tot,
        "u_L": u_L,
        "d2F_sym": d2F_sym,
    }


def symbolic():
    import sympy as sp
    W, V_ratio, f_V, gamma, a = sp.symbols("W V_ratio f_V gamma a", positive=True)
    x1 = (f_V * V_ratio) ** sp.Rational(1, 3)
    x2 = ((1 - f_V) * V_ratio) ** sp.Rational(1, 3)
    F1 = x1 ** 2 - W * (1 - 1 / sp.sqrt(1 + x1 ** 2))
    F2 = x2 ** 2 - W * (1 - 1 / sp.sqrt(1 + x2 ** 2))
    F_tot_dimless = F1 + F2
    x_s = (V_ratio / 2) ** sp.Rational(1, 3)
    d2F_sym = -4 + 2 * W * (1 + 4 * x_s ** 2) / (1 + x_s ** 2) ** sp.Rational(5, 2)
    return {
        "F_tot_dimless": F_tot_dimless,
        "F_tot": (gamma * sp.Rational(1, 1000)) * a ** 2 * F_tot_dimless,
        "u_L": x1 / (x1 + x2),
        "d2F_sym": d2F_sym,
    }
