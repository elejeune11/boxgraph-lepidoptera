import math

K_B = 1.380649e-23  # J/K (exact SI)


def run(inputs, settings):
    R_m = inputs["R"] * 1.0e-6
    lp_m = inputs["l_p"] * 1.0e-6
    d_m = inputs["d"] * 1.0e-9
    gamma_si = inputs["gamma"] * 1.0e-3
    T = inputs["T"]
    alpha = (K_B * T * lp_m) / (2.0 * math.pi * d_m * gamma_si * (R_m ** 2))
    L_star_m = (2.0 * R_m) / (1.0 - alpha)
    return {"L_star": L_star_m * 1.0e6}


def symbolic():
    import sympy as sp

    R, l_p, d, gamma, T = sp.symbols("R l_p d gamma T", positive=True)
    alpha = sp.Float(K_B * 1.0e18) * T * l_p / (2 * sp.pi * d * gamma * R ** 2)
    return {"L_star": (2 * R) / (1 - alpha)}
