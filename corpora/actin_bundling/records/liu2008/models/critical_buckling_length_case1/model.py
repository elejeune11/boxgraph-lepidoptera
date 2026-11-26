import math


def run(inputs, settings):
    kBT = float(inputs["kBT"])      # pN*nm
    L_p = float(inputs["L_p"])      # um
    f_mem = float(inputs["f_mem"])  # pN
    # kBT / f_mem is in nm -> * 1e-3 for um; product with L_p (um) is in um^2
    L_buckle = math.sqrt((math.pi ** 2) * (kBT * 1e-3) * L_p / (4.0 * f_mem))
    return {"L_buckle": L_buckle}


def symbolic():
    import sympy as sp

    kBT, L_p, f_mem = sp.symbols("kBT L_p f_mem", positive=True)
    return {"L_buckle": sp.sqrt(sp.pi**2 * kBT * sp.Rational(1, 1000) * L_p / (4 * f_mem))}
