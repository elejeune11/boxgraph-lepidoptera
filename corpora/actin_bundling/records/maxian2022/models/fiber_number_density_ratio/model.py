def run(inputs, settings):
    rho_short = inputs["F_short"] / (inputs["Ld_short"] ** 3)
    rho_long = inputs["F_long"] / (inputs["Ld_long"] ** 3)
    return {
        "rho_short": rho_short,
        "rho_long": rho_long,
        "ratio": rho_short / rho_long,
    }


def symbolic():
    import sympy as sp

    F_short, Ld_short, F_long, Ld_long = sp.symbols("F_short Ld_short F_long Ld_long", positive=True)
    rho_short = F_short / Ld_short**3
    rho_long = F_long / Ld_long**3
    return {
        "rho_short": rho_short,
        "rho_long": rho_long,
        "ratio": rho_short / rho_long,
    }
