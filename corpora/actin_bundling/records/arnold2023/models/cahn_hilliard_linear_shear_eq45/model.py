import math
import numpy as np


def run(inputs, settings):
    phi0 = float(inputs["phi_0"])
    Pe = float(inputs["Pe"])
    t = float(inputs["t"])
    n_k = int(settings.get("n_k", 2000))
    d_ln_t = float(settings.get("d_ln_t", 0.01))

    c = 1.0 - 3.0 * phi0 * phi0

    def compute_a(t_val, Pe_val, ray_sign):
        k_scale = t_val ** (-0.25)
        k_hi = max(6.0 * k_scale, math.sqrt(max(c, 0.0)) + 6.0 * k_scale)
        ks = np.linspace(1.0e-5 * k_scale, k_hi, n_k)
        wt = (c * ks * ks - ks ** 4) * t_val
        e2w = np.exp(2.0 * (wt - np.max(wt)))
        if ray_sign == 0.0:
            # 2D angular average of S(k, theta, t) to O(Pe): <k_x k_y>_theta = (k^2/2)<sin(2*theta)> = 0
            S = e2w
        else:
            # Single diagonal ray k_x = ray_sign * k_y (k_x k_y = ray_sign * k^2 / 2) to O(Pe)
            kxky = ray_sign * 0.5 * ks * ks
            S = e2w * (1.0 + 2.0 * Pe_val * kxky * (c - 2.0 * ks * ks) * (t_val ** 2))
        return float(np.trapezoid(S, ks) / np.trapezoid(ks * S, ks))

    t2 = t * math.exp(d_ln_t)
    a_rad_1 = compute_a(t, Pe, 0.0)
    a_rad_0 = compute_a(t, 0.0, 0.0)
    a_rad_2 = compute_a(t2, 0.0, 0.0)
    alpha_passive = math.log(a_rad_2 / a_rad_0) / d_ln_t

    # O(Pe) perturbation magnitude along diagonal ray k_x = +k_y (ray_sign = +1 gives da_ray > 0 for c=0)
    da_ray_1 = compute_a(t, Pe, 1.0) - compute_a(t, 0.0, 1.0)
    da_ray_2 = compute_a(t2, Pe, 1.0) - compute_a(t2, 0.0, 1.0)
    alpha_shear_ray_correction = math.log(abs(da_ray_2) / abs(da_ray_1)) / d_ln_t

    return {
        "a_rad": float(a_rad_1),
        "alpha_passive": float(alpha_passive),
        "alpha_shear_ray_correction": float(alpha_shear_ray_correction),
        "delta_a_rad_O_Pe": float(a_rad_1 - a_rad_0),
    }
