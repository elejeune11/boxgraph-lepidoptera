import math


def run(inputs, settings):
    chi = float(inputs["chi_mu"])
    phi_g = 0.0 if chi >= 0.0 else math.sqrt(-3.0 * chi)
    return {"phi_g": phi_g}


def symbolic():
    import sympy as sp

    chi_mu = sp.Symbol("chi_mu", real=True)
    return {"phi_g": sp.Piecewise((0, chi_mu >= 0), (sp.sqrt(-3 * chi_mu), True))}
