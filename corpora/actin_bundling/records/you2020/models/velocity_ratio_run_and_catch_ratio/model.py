def run(inputs, settings):
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta_over_kappa"]) * k
    return {"v_ratio_fs": (d * d - k * k) / (aB * aB)}


def symbolic():
    import sympy as sp

    alpha_B, kappa, delta_over_kappa = sp.symbols("alpha_B kappa delta_over_kappa", real=True)
    d = delta_over_kappa * kappa
    return {"v_ratio_fs": (d**2 - kappa**2) / alpha_B**2}
