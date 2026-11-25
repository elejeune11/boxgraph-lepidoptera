import math


def run(inputs, settings):
    K = inputs["K"]  # pN
    gamma = inputs["gamma"] * 1e-3  # nN/m -> pN/um
    R = inputs["R"]  # um
    a = inputs["alpha"]  # rad

    sa = math.sin(a)
    ca = math.cos(a)
    f_el = -math.pi * (7.0 * a * ca - 7.0 * sa + (a ** 2) * sa)
    f_s = 4.0 * math.pi * (sa - a * ca)
    v = 2.0 * math.pi * (sa - a * ca - (sa ** 3) / 3.0)

    V = (R ** 3) * v
    F_el = K * R * f_el
    F_s = gamma * (R ** 2) * f_s
    return {
        "f_el": f_el,
        "f_s": f_s,
        "v": v,
        "V": V,
        "F_el": F_el,
        "F_s": F_s,
        "F_bare": F_el + F_s,
    }


def symbolic():
    import sympy as sp
    K, gamma, R, alpha = sp.symbols("K gamma R alpha", positive=True)
    gamma_pN = gamma * sp.Rational(1, 1000)
    f_el = -sp.pi * (7 * alpha * sp.cos(alpha) - 7 * sp.sin(alpha) + alpha ** 2 * sp.sin(alpha))
    f_s = 4 * sp.pi * (sp.sin(alpha) - alpha * sp.cos(alpha))
    v = 2 * sp.pi * (sp.sin(alpha) - alpha * sp.cos(alpha) - sp.sin(alpha) ** 3 / 3)
    F_el = K * R * f_el
    F_s = gamma_pN * R ** 2 * f_s
    return {
        "f_el": f_el,
        "f_s": f_s,
        "v": v,
        "V": R ** 3 * v,
        "F_el": F_el,
        "F_s": F_s,
        "F_bare": F_el + F_s,
    }
