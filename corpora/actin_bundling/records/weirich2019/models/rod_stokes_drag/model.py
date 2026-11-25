import math


def run(inputs, settings):
    # 1 Pa*s = 1 N*s/m^2 = 1e12 pN*s / (1e12 um^2) = 1 pN*s/um^2
    eta = inputs["eta"]              # pN*s/um^2
    L_M = inputs["L_M"]              # um
    D_M_um = inputs["D_M"] * 1e-3    # nm -> um
    Gamma = 6.0 * math.pi * eta * L_M / math.log(L_M / D_M_um)
    return {"Gamma": Gamma}


def symbolic():
    import sympy as sp
    eta, L_M, D_M = sp.symbols("eta L_M D_M", positive=True)
    return {"Gamma": 6 * sp.pi * eta * L_M / sp.log(L_M / (D_M * sp.Rational(1, 1000)))}
