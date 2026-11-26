import math


def run(inputs, settings):
    n_nat = float(inputs["n_native"])
    t_nat = float(inputs["t_native"])
    n_bun = float(inputs["n_bundle"])
    t_bun = float(inputs["t_bundle"])
    omega_0 = 2.0 * math.pi * t_nat / n_nat
    omega_bundle = 2.0 * math.pi * t_bun / n_bun
    delta_omega_deg = (omega_bundle - omega_0) * (180.0 / math.pi)
    return {
        "omega_0": omega_0,
        "omega_bundle": omega_bundle,
        "delta_omega": delta_omega_deg,
    }


def symbolic():
    import sympy as sp

    n_native, t_native, n_bundle, t_bundle = sp.symbols(
        "n_native t_native n_bundle t_bundle", positive=True
    )
    omega_0 = 2 * sp.pi * t_native / n_native
    omega_bundle = 2 * sp.pi * t_bundle / n_bundle
    return {
        "omega_0": omega_0,
        "omega_bundle": omega_bundle,
        "delta_omega": (omega_bundle - omega_0) * 180 / sp.pi,
    }
