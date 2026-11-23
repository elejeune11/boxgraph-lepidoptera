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


def run(inputs, settings):
    l_dom = float(inputs["l"])
    xi = float(inputs["xi"])
    chi = float(inputs["chi_rho0"])
    eps = float(inputs["epsilon"])
    noise_amp = float(inputs["noise_amp"])
    t_total = float(inputs["t_total"])

    N_grid = int(round(settings.get("N_grid", 512)))
    dt = float(settings.get("dt", 0.05))
    seed = int(round(settings.get("seed", 1)))

    rng = np.random.default_rng(seed)
    dx = l_dom / N_grid
    k = np.arange(N_grid)
    lam_lap = (4.0 / (dx * dx)) * (np.sin(np.pi * k / (2.0 * N_grid)) ** 2)
    helm_inv = 1.0 / (1.0 + lam_lap / xi)

    rho = np.maximum(1.0 + noise_amp * rng.standard_normal(N_grid), 1e-6)
    n_steps = max(1, int(round(t_total / dt)))
    sample_stride = max(1, int(round(1.0 / dt)))  # sample at unit time intervals dt_sample = 1.0
    c_fac = dt / (dx * dx)
    ab = np.zeros((3, N_grid))

    peak_counts = []
    min_rho = float(np.min(rho))
    max_rho = float(np.max(rho))

    for step in range(n_steps):
        s = dct(dct(chi * (rho - 1.0), type=2, norm="ortho") * helm_inv, type=3, norm="ortho")
        ds_pe = (s[1:] - s[:-1]) / xi
        bp = c_fac * _bernoulli(ds_pe)
        bm = c_fac * _bernoulli(-ds_pe)

        diag = np.full(N_grid, 1.0 + dt)
        upper = np.zeros(N_grid)
        lower = np.zeros(N_grid)
        diag[:-1] += bm
        upper[1:] -= bp
        diag[1:] += bp
        lower[:-1] -= bm

        ab[0, :] = upper
        ab[1, :] = diag
        ab[2, :] = lower

        # Semi-implicit treatment of -eps*(rho - 1)^2:
        # keep positivity by linearizing quadratic decay when rho > 1
        rhs = rho + dt * (1.0 - eps * (rho - 1.0) ** 2)
        rho = solve_banded((1, 1), ab, rhs, overwrite_ab=False, overwrite_b=True, check_finite=False)
        mn = float(np.min(rho))
        mx = float(np.max(rho))
        if mn < min_rho:
            min_rho = mn
        if mx > max_rho:
            max_rho = mx

        if step % sample_stride == 0:
            thresh = mn + 0.1 * (mx - mn)
            above = rho > thresh
            pks = int(np.sum(above & np.concatenate(([True], ~above[:-1]))))
            peak_counts.append(pks)

    second_half = peak_counts[len(peak_counts) // 2 :]
    diffs = np.diff(second_half) if len(second_half) > 1 else np.array([0])
    n_events = int(np.sum(diffs != 0))
    return {
        "n_events": float(n_events),
        "min_rho": min_rho,
        "max_rho": max_rho,
    }
