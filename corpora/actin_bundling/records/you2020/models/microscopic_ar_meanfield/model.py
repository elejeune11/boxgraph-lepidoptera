import math


def run(inputs, settings):
    kr = float(inputs["k_r"])
    sr = float(inputs["sigma_r"])
    kmn = float(inputs["k_munu"])
    smn = float(inputs["sigma_munu"])
    phi = float(inputs["phi_mu"])
    pref = math.pi / 6.0
    R = (kr * sr**4 + kmn * smn**4) * pref
    kappa_eff = kr * sr**4 * pref * phi
    delta_eff = kmn * smn**4 * pref * phi
    return {"R_munu": R, "kappa_eff": kappa_eff, "delta_eff": delta_eff}


def symbolic():
    import sympy as sp

    k_r, sigma_r, k_munu, sigma_munu, phi_mu = sp.symbols(
        "k_r sigma_r k_munu sigma_munu phi_mu", real=True
    )
    pref = sp.pi / 6
    return {
        "R_munu": (k_r * sigma_r**4 + k_munu * sigma_munu**4) * pref,
        "kappa_eff": k_r * sigma_r**4 * pref * phi_mu,
        "delta_eff": k_munu * sigma_munu**4 * pref * phi_mu,
    }
