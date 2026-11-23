import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import root


def run(inputs, settings):
    alpha = float(inputs["alpha"])
    phi = float(inputs["phi"])
    i_mode = int(round(inputs["i"]))
    x_targets = [float(v) for v in inputs["target_branch"]["x"]]
    x_max_need = max(x_targets) * 1.03

    M = int(round(settings.get("M_modes", 40)))
    N = int(round(settings.get("N_quad", 512)))
    n_steps = int(round(settings.get("n_steps", 70)))
    if N < 8 * M:
        N = 8 * M

    # Single-period phase coordinate theta = i_mode * u in [0, pi] gives full
    # spectral accuracy per wavelength regardless of mode number i_mode.
    j = np.arange(N)
    theta = np.pi * (j + 0.5) / N
    m_idx = np.arange(1, M + 1)
    harms = m_idx * i_mode
    C = np.cos(np.outer(theta, m_idx))
    S = np.sin(np.outer(theta, m_idx))
    alt_sign = (-1.0) ** m_idx
    k_n = harms * np.pi / alpha

    k_i = i_mode * np.pi
    chi0 = 1.0 + (k_i / alpha) ** 2 + phi * (1.0 + (alpha / k_i) ** 2)
    z = np.zeros(M)
    z[0] = chi0

    # Continue monotonically in first-harmonic amplitude eps until max(rho) >= x_max_need
    max_rhos = [1.0]
    chis = [chi0]
    min_rhos = [1.0]
    eps = 0.0
    deps = 1.45 / n_steps
    for _ in range(n_steps * 2):
        eps += deps

        def resid(zv):
            chi_v = zv[0]
            r_c = np.zeros(M)
            r_c[0] = eps
            r_c[1:] = zv[1:]
            s_c = chi_v * r_c / (1.0 + k_n * k_n)
            rho_g = 1.0 + C @ r_c
            ds_g = -(S @ (harms * np.pi * s_c))
            flux_s = (2.0 / N) * (S.T @ (rho_g * ds_g))
            return (harms * np.pi / (alpha * alpha)) * flux_s + (k_n * k_n + phi) * r_c

        sol = root(resid, z, tol=1e-12)
        if not sol.success:
            eps -= deps
            deps *= 0.5
            if deps < 1e-4:
                break
            continue
        z = sol.x
        chi_v = float(z[0])
        r_c = np.zeros(M)
        r_c[0] = eps
        r_c[1:] = z[1:]
        rho_g = 1.0 + C @ r_c
        rho0 = float(1.0 + np.sum(r_c))
        rhopi = float(1.0 + np.sum(r_c * alt_sign))
        mx = max(rho0, float(np.max(rho_g)))
        mn = min(rhopi, float(np.min(rho_g)))
        max_rhos.append(mx)
        chis.append(chi_v)
        min_rhos.append(mn)
        if mx >= x_max_need:
            break

    max_rhos_arr = np.array(max_rhos)
    chis_arr = np.array(chis)
    min_rhos_arr = np.array(min_rhos)
    cs_chi = CubicSpline(max_rhos_arr, chis_arr)
    cs_min = CubicSpline(max_rhos_arr, min_rhos_arr)
    y_preds = [float(cs_chi(xv)) for xv in x_targets]
    x_hi = max(x_targets)
    return {
        "chi_curve": {"x": x_targets, "y": y_preds},
        "chi_at_max": float(cs_chi(x_hi)),
        "min_rho_at_max": float(cs_min(x_hi)),
    }
