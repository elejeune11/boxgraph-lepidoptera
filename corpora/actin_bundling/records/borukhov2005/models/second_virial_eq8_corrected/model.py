import math

_EULER_GAMMA = 0.5772156649015329


def run(inputs, settings):
    gamma = inputs["gamma"]
    D = inputs["D"]
    L = inputs["L"]
    delta = inputs["delta"]
    kappa_inv = inputs["kappa_inv"]
    Gamma_par = inputs["Gamma_par"]
    Gamma_perp = inputs["Gamma_perp"]
    u_par = inputs["u_par"]
    u_perp = inputs["u_perp"]

    D_eff = D + kappa_inv * (_EULER_GAMMA + math.log(Gamma_perp))

    if gamma <= D / L:
        t1 = 8.0 * (D ** 2 + (kappa_inv ** 2) * (math.log(Gamma_par * L / D) ** 2)) / (D_eff * L)
        Y = (Gamma_par - abs(u_par)) * (L / D)
        if abs(Y) < 1e-8:
            bracket = 0.5 * Y
        else:
            bracket = 1.0 + (math.exp(-Y) - 1.0) / Y
        t2 = 16.0 * (D * delta / (D_eff * L)) * bracket
        return {"B": float(t1 + t2)}
    else:
        sg = math.sin(gamma)
        t1 = (8.0 / math.pi) * sg * (1.0 - (kappa_inv / D_eff) * math.log(sg))
        t2 = (8.0 * delta / (math.pi * D_eff)) * sg * (1.0 - math.exp(-Gamma_perp / sg + abs(u_perp)))
        return {"B": float(t1 + t2)}
