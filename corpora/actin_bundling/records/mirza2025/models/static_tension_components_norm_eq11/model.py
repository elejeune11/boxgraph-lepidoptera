def run(inputs, settings):
    kap = float(inputs["kappa"])
    S = float(inputs["S"])
    _ = float(inputs["t_norm"])  # static steady-state prediction is constant in time
    return {
        "sigma_par_norm": 1.0 + 0.5 * kap * S,
        "sigma_perp_norm": 1.0 - 0.5 * kap * S,
    }


def symbolic():
    import sympy as sp
    kappa, S, t_norm = sp.symbols("kappa S t_norm", real=True)
    return {
        "sigma_par_norm": 1 + kappa * S / 2,
        "sigma_perp_norm": 1 - kappa * S / 2,
    }
