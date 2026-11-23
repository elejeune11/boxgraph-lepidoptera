import numpy as np
from scipy.optimize import root


def run(inputs, settings):
    alpha = float(inputs["alpha"])
    phi = float(inputs["phi"])
    i_mode = int(round(inputs["i"]))
    eps = float(inputs["eps"])
    M = int(round(settings.get("M_modes", 24)))
    N = int(round(settings.get("N_quad", 256)))

    if M <= 2 * i_mode:
        M = 2 * i_mode + 4
    if N < 4 * M:
        N = 4 * M

    j = np.arange(N)
    u_grid = np.pi * (j + 0.5) / N
    modes = np.arange(M)
    C = np.cos(np.outer(u_grid, modes))
    S = np.sin(np.outer(u_grid, modes))
    k_n = modes * np.pi / alpha

    idx_other = [n for n in range(1, M) if n != i_mode]

    def residual(z):
        chi_v = z[0]
        r_c = np.zeros(M)
        r_c[i_mode] = eps
        for idx, n in enumerate(idx_other):
            r_c[n] = z[1 + idx]
        s_c = chi_v * r_c / (1.0 + k_n * k_n)
        rho_g = 1.0 + C @ r_c
        ds_g = -(S @ (modes * np.pi * s_c))
        flux_s = (2.0 / N) * (S.T @ (rho_g * ds_g))
        eq = np.empty(M - 1)
        for n in range(1, M):
            eq[n - 1] = (n * np.pi / (alpha * alpha)) * flux_s[n] + (k_n[n] ** 2 + phi) * r_c[n]
        return eq

    k_i = i_mode * np.pi
    chi0 = 1.0 + (k_i / alpha) ** 2 + phi * (1.0 + (alpha / k_i) ** 2)
    z0 = np.zeros(M - 1)
    z0[0] = chi0
    sol = root(residual, z0, tol=1e-12)
    if not sol.success:
        raise RuntimeError(f"steady-state bifurcation solver failed: {sol.message}")

    chi_v = float(sol.x[0])
    r_c = np.zeros(M)
    r_c[i_mode] = eps
    for idx, n in enumerate(idx_other):
        r_c[n] = sol.x[1 + idx]
    rho_g = 1.0 + C @ r_c

    return {
        "chi_branch": chi_v,
        "chi_i2": 2.0 * (chi_v - chi0) / (eps * eps),
        "min_rho": float(np.min(rho_g)),
    }
