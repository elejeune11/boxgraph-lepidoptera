import numpy as np
from scipy.fft import dct
from scipy.linalg import solve_banded


def _bernoulli(z):
    out = np.empty_like(z, dtype=float)
    small = np.abs(z) < 1e-4
    zs = z[small]
    out[small] = 1.0 - 0.5 * zs + (zs * zs) / 12.0 - (zs ** 4) / 720.0
    zl = z[~small]
    pos = zl > 500.0
    neg = zl < -500.0
    mid = ~(pos | neg)
    res = np.empty_like(zl)
    res[pos] = 0.0
    res[neg] = -zl[neg]
    zm = zl[mid]
    res[mid] = zm / np.expm1(zm)
    out[~small] = res
    return out


def _single_realization(alpha, phi, chi_start, chi_end, t_ramp_nd, t_tot_nd, N_grid, dt, noise_amp, rng):
    dy = 1.0 / N_grid
    k = np.arange(N_grid)
    lam_lap = (4.0 / (dy * dy)) * (np.sin(np.pi * k / (2.0 * N_grid)) ** 2)
    helm_inv = 1.0 / (1.0 + lam_lap / (alpha * alpha))

    rho = np.maximum(1.0 + noise_amp * rng.standard_normal(N_grid), 1e-6)
    n_steps = max(1, int(round(t_tot_nd / dt)))
    c_fac = dt / (alpha * alpha * dy * dy)

    ab = np.zeros((3, N_grid))
    min_rho = float(np.min(rho))

    for step in range(n_steps):
        t = step * dt
        if t_ramp_nd > 0.0 and t < t_ramp_nd:
            chi = chi_start + (chi_end - chi_start) * (t / t_ramp_nd)
        else:
            chi = chi_end

        s = dct(dct(chi * (rho - 1.0), type=2, norm="ortho") * helm_inv, type=3, norm="ortho")
        ds = s[1:] - s[:-1]
        bp = c_fac * _bernoulli(ds)
        bm = c_fac * _bernoulli(-ds)

        diag = np.full(N_grid, 1.0 + dt * phi)
        upper = np.zeros(N_grid)
        lower = np.zeros(N_grid)
        diag[:-1] += bm
        upper[1:] -= bp
        diag[1:] += bp
        lower[:-1] -= bm

        ab[0, :] = upper
        ab[1, :] = diag
        ab[2, :] = lower

        rho = solve_banded((1, 1), ab, rho + dt * phi, overwrite_ab=False, overwrite_b=True, check_finite=False)
        mn = float(np.min(rho))
        if mn < min_rho:
            min_rho = mn

    S = float(np.mean((rho - 1.0) ** 2))
    # Peak definition from SI Section 2.2 (p. 11): > 10% of (max - min) above min
    thresh = np.min(rho) + 0.1 * (np.max(rho) - np.min(rho))
    above = rho > thresh
    peaks = int(np.sum(above & np.concatenate(([True], ~above[:-1]))))
    return float(peaks), S, min_rho


def run(inputs, settings):
    alpha = float(inputs["alpha"])
    phi = float(inputs["phi"])
    tau = float(inputs["tau"])
    chi_start = float(inputs["chi_start"])
    chi_end = float(inputs["chi_end"])
    t_ramp = float(inputs["t_ramp"])
    t_total = float(inputs["t_total"])
    noise_amp = float(inputs["noise_amp"])

    N_grid = int(round(settings.get("N_grid", 2048)))
    dt = float(settings.get("dt", 100.0))
    n_samples = max(1, int(round(settings.get("n_samples", 4))))
    seed = int(round(settings.get("seed", 0)))

    # Nondimensional time scale T = eta / (xi * D) = phi * tau
    T_scale = phi * tau
    t_ramp_nd = t_ramp / T_scale
    t_tot_nd = t_total / T_scale

    rng = np.random.default_rng(seed)
    peaks_list = []
    S_list = []
    min_rhos = []
    for _ in range(n_samples):
        pk, S, mn = _single_realization(
            alpha, phi, chi_start, chi_end, t_ramp_nd, t_tot_nd, N_grid, dt, noise_amp, rng
        )
        peaks_list.append(pk)
        S_list.append(S)
        min_rhos.append(mn)

    return {
        "n_rings": float(np.mean(peaks_list)),
        "S": float(np.mean(S_list)),
        "min_rho": float(np.min(min_rhos)),
    }
