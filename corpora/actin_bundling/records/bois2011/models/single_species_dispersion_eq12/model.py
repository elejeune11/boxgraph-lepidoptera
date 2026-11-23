def run(inputs, settings):
    k = inputs["k"]
    D = inputs["D"]
    Pe = inputs["Pe"]
    c0 = inputs["c0"]
    df_dc = inputs["df_dc"]
    l = inputs["l"]
    kl2 = (k * l) ** 2
    factor = 1.0 - Pe * c0 * df_dc / (1.0 + kl2)
    lam = -(k**2) * D * factor
    lam_dimless = -kl2 * factor
    return {"lam": lam, "lam_dimless": lam_dimless}


def symbolic():
    import sympy as sp

    k, D, Pe, c0, df_dc, l = sp.symbols("k D Pe c0 df_dc l", positive=True)
    lam = -(k**2) * D * (1 - Pe * c0 * df_dc / (1 + k**2 * l**2))
    lam_dimless = -(k**2 * l**2) * (1 - Pe * c0 * df_dc / (1 + k**2 * l**2))
    return {"lam": lam, "lam_dimless": lam_dimless}
