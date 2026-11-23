from scipy.optimize import brentq


def run(inputs, settings):
    R = inputs["R"]
    mu_par_bar = inputs["mu_par_bar"]
    mu_perp_bar = inputs["mu_perp_bar"]
    K3_bar = inputs["K3_bar"]
    gamma = inputs["gamma"]
    rtol = float(settings.get("rtol", 1e-8))

    a = 4.0 * (mu_perp_bar * (R**4) + K3_bar * (R**2))
    b = 2.0 * mu_par_bar * (R**2)
    hi = min(gamma / b, (gamma / a) ** (1.0 / 3.0)) if b > 0.0 else (gamma / a) ** (1.0 / 3.0)
    om0 = brentq(lambda om: a * (om**3) + b * om - gamma, 0.0, hi * (1.0 + 1e-12), rtol=rtol, xtol=1e-300)
    e_elast = (
        mu_par_bar * ((om0 * R) ** 2)
        + mu_perp_bar * ((om0 * R) ** 4)
        + K3_bar * (om0**4) * (R**2)
        - gamma * om0
    )
    return {"Omega_0": om0, "theta_0": om0 * R, "E_elast": e_elast}
