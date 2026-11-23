def run(inputs, settings):
    k2 = inputs["k"] ** 2
    D = inputs["D"]
    kl2 = k2 * (inputs["l"] ** 2)
    lam = (
        -k2 * D
        + inputs["dR_dc"]
        + D * k2 * inputs["Pe"] * inputs["c0"] * inputs["df_dc"] / (1.0 + kl2)
    )
    return {"lam": lam}


def symbolic():
    import sympy as sp

    k, D, Pe, c0, df_dc, l = sp.symbols("k D Pe c0 df_dc l", positive=True)
    dR_dc = sp.Symbol("dR_dc", real=True)
    return {"lam": -(k**2) * D + dR_dc + D * k**2 * Pe * c0 * df_dc / (1 + k**2 * l**2)}
