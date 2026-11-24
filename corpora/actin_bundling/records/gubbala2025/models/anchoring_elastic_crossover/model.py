import math


def run(inputs, settings):
    chi = float(inputs["chi"])
    R = float(inputs["R"])
    E1 = float(inputs["E1"])
    E3 = float(inputs["E3"])
    S = float(inputs["S"])

    # 1D bare Cahn-Hilliard kink integral: int (d_s phi_0)^2 ds = 2*sqrt(2)/3,
    # so f_anch = (chi/2)*(n . grad phi)^2 integrates to gamma_A_tilde = (sqrt(2)/3)*chi.
    gamma_A_tilde = (math.sqrt(2.0) / 3.0) * chi
    num = gamma_A_tilde * R

    k11_eq6 = 2.0 * E1 * S**2 - (2.0 / 3.0) * E3 * S**3
    k33_eq6 = 2.0 * E1 * S**2 + (4.0 / 3.0) * E3 * S**3
    k11_2d = 2.0 * E1 * S**2 - E3 * S**3
    k33_2d = 2.0 * E1 * S**2 + E3 * S**3

    return {
        "xi_splay_eq6": num / k11_eq6,
        "xi_bend_eq6": num / k33_eq6,
        "xi_splay_2d": num / k11_2d,
        "xi_bend_2d": num / k33_2d,
    }


def symbolic():
    import sympy as sp

    chi, R, E1, E3, S = sp.symbols("chi R E1 E3 S", positive=True)
    num = sp.sqrt(2) / 3 * chi * R
    k11_eq6 = 2 * E1 * S**2 - sp.Rational(2, 3) * E3 * S**3
    k33_eq6 = 2 * E1 * S**2 + sp.Rational(4, 3) * E3 * S**3
    k11_2d = 2 * E1 * S**2 - E3 * S**3
    k33_2d = 2 * E1 * S**2 + E3 * S**3
    return {
        "xi_splay_eq6": num / k11_eq6,
        "xi_bend_eq6": num / k33_eq6,
        "xi_splay_2d": num / k11_2d,
        "xi_bend_2d": num / k33_2d,
    }
