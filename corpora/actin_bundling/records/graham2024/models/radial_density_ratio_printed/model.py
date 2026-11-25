def run(inputs, settings):
    num = inputs["N_shell"] / (inputs["r"] ** 2 * inputs["delta_r"])
    den = inputs["N_total"] / (inputs["R_drop"] ** 3)
    return {"density_ratio": num / den}


def symbolic():
    import sympy as sp

    N_shell, r, delta_r, N_total, R_drop = sp.symbols(
        "N_shell r delta_r N_total R_drop", positive=True
    )
    return {"density_ratio": (N_shell / (r**2 * delta_r)) / (N_total / R_drop**3)}
