import math


def run(inputs, settings):
    r = inputs["r"]
    gamma = inputs["gamma"]
    D = inputs["D"]
    L = inputs["L"]
    delta = inputs["delta"]
    kappa_inv = inputs["kappa_inv"]
    u_par = inputs["u_par"]
    u_perp = inputs["u_perp"]
    Gamma_par = inputs["Gamma_par"]
    Gamma_perp = inputs["Gamma_perp"]

    if r <= D:
        return {"u": float("inf")}

    if gamma <= D / L:
        u_g = (L / D) * u_par
        G_g = (L / D) * Gamma_par
    else:
        sg = math.sin(gamma)
        u_g = u_perp
        G_g = Gamma_perp / sg

    if r <= D + delta:
        val = u_g + G_g
    else:
        val = G_g * math.exp(-(r - D - delta) / kappa_inv)
    return {"u": float(val)}
