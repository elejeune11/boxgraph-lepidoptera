import math

KB_PN_UM_PER_K = 1.380649e-5


def run(inputs, settings):
    # 1 Pa*s = 1 pN*s/um^2, so 1/(6*pi*R[um]*nu[Pa*s]) is directly in um/(pN*s).
    mu = 1.0 / (6.0 * math.pi * inputs["R"] * inputs["nu"])
    D = KB_PN_UM_PER_K * inputs["T"] * mu
    return {"mu": mu, "D": D}


def symbolic():
    import sympy as sp

    R, nu, T = sp.symbols("R nu T", positive=True)
    mu = 1 / (6 * sp.pi * R * nu)
    return {"mu": mu, "D": KB_PN_UM_PER_K * T * mu}
