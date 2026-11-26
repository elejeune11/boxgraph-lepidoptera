import math


def run(inputs, settings):
    kBT = float(inputs["kBT"])  # pN*nm
    L_p = float(inputs["L_p"])  # um
    L = float(inputs["L"])      # um
    # (pN*nm) * um / um^2 = pN * (nm / um) = 1e-3 pN = 1 fN
    f_buckle = (math.pi ** 2) * kBT * L_p / (4.0 * L ** 2)
    return {"f_buckle": f_buckle}


def symbolic():
    import sympy as sp

    kBT, L_p, L = sp.symbols("kBT L_p L", positive=True)
    return {"f_buckle": sp.pi**2 * kBT * L_p / (4 * L**2)}
