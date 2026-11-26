import math

_EULER_GAMMA = 0.5772156649015329


def run(inputs, settings):
    D = inputs["D"]
    kappa_inv = inputs["kappa_inv"]
    Gamma_perp = inputs["Gamma_perp"]
    D_eff = D + kappa_inv * (_EULER_GAMMA + math.log(Gamma_perp))
    return {"D_eff": float(D_eff)}


def symbolic():
    import sympy as sp

    D, kappa_inv, Gamma_perp = sp.symbols("D kappa_inv Gamma_perp", positive=True)
    return {"D_eff": D + kappa_inv * (sp.EulerGamma + sp.log(Gamma_perp))}
