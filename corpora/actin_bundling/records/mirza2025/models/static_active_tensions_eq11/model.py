def run(inputs, settings):
    rho = inputs["rho"]
    rho_ref = inputs["rho_ref"]
    lam = inputs["lam"]
    kap = inputs["kappa"]
    S = inputs["S"]
    sigma_par = rho * lam * (1.0 + 0.5 * kap * S)
    sigma_perp = rho * lam * (1.0 - 0.5 * kap * S)
    sigma_bar = 0.5 * (sigma_par + sigma_perp)
    sigma_bar_ref = rho_ref * lam
    sigma_bar_norm = sigma_bar / sigma_bar_ref
    rho_norm = rho / rho_ref
    norm_slope = sigma_bar_norm / rho_norm
    xi = (sigma_par - sigma_perp) / sigma_bar
    return {
        "sigma_par": sigma_par,
        "sigma_perp": sigma_perp,
        "sigma_bar": sigma_bar,
        "sigma_bar_norm": sigma_bar_norm,
        "norm_slope": norm_slope,
        "xi": xi,
    }

def symbolic():
    import sympy as sp
    rho, rho_ref, lam, kappa, S = sp.symbols("rho rho_ref lam kappa S", real=True)
    sigma_par = rho * lam * (1 + kappa * S / 2)
    sigma_perp = rho * lam * (1 - kappa * S / 2)
    sigma_bar = (sigma_par + sigma_perp) / 2
    sigma_bar_norm = sigma_bar / (rho_ref * lam)
    norm_slope = sp.simplify(sigma_bar_norm / (rho / rho_ref))
    xi = sp.simplify((sigma_par - sigma_perp) / sigma_bar)
    return {
        "sigma_par": sigma_par,
        "sigma_perp": sigma_perp,
        "sigma_bar": sigma_bar,
        "sigma_bar_norm": sigma_bar_norm,
        "norm_slope": norm_slope,
        "xi": xi,
    }
