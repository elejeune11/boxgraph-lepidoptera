import math


def run(inputs, settings):
    geom = 1.0 + (inputs["beta"] * math.pi / inputs["L_over_l"]) ** 2
    denom = inputs["c0"] * inputs["df_dc"]
    return {
        "stability_ratio": inputs["Pe"] * denom / geom,
        "Pe_crit": geom / denom,
    }


def symbolic():
    import sympy as sp

    Pe, c0, df_dc, beta, L_over_l = sp.symbols("Pe c0 df_dc beta L_over_l", positive=True)
    geom = 1 + (beta * sp.pi / L_over_l) ** 2
    return {
        "stability_ratio": Pe * c0 * df_dc / geom,
        "Pe_crit": geom / (c0 * df_dc),
    }
