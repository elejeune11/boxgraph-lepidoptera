import math
import numpy as np
from scipy.optimize import root


def _solve_steady_n1_profile(a, a_Delta, kappa, rho, alpha, ell, p1_bar, p2_bar, Nx):
    M = min(16, Nx // 4)
    N = 4 * M
    theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    ms = np.arange(1, M + 1)
    C = np.cos(np.outer(theta, ms))
    k_sq = (2.0 * np.pi * ms / ell) ** 2

    def res(u):
        c1, c2 = u[:M], u[M:]
        p1 = p1_bar + C @ c1
        p2 = p2_bar + C @ c2
        f1p = a * p1 + p1**3 - (rho + alpha) * p2
        f2p = (a + a_Delta) * p2 + p2**3 - (rho - alpha) * p1
        r1 = k_sq * c1 + (2.0 / N) * (C.T @ f1p)
        r2 = kappa * k_sq * c2 + (2.0 / N) * (C.T @ f2p)
        return np.concatenate([r1, r2])

    u0 = np.zeros(2 * M)
    u0[0] = 0.3
    u0[M] = 0.3
    sol = root(res, u0, method="hybr", tol=1e-12)
    c1, c2 = sol.x[:M], sol.x[M:]
    x = np.linspace(0.0, 1.0, Nx, endpoint=False)
    Cx = np.cos(np.outer(2.0 * np.pi * x, ms))
    return p1_bar + Cx @ c1, p2_bar + Cx @ c2


def run(inputs, settings):
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    Q = inputs["Q"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    ell = inputs["ell"]
    p1_bar = inputs["phi1_bar"]
    p2_bar = inputs["phi2_bar"]
    noise_amp = inputs["noise_amp"]
    t_end = inputs["t_end"]

    Nx = int(settings.get("n_grid", 64))
    dt = float(settings.get("dt", 0.1))
    S = float(settings.get("S_stab", 6.0))
    seed = int(settings.get("seed", 0))
    init_type = str(settings.get("init_type", "homogeneous"))

    rng = np.random.default_rng(seed)
    # Generate white noise in low-frequency Fourier modes so IC is grid-consistent under Nx refinement
    n_noise_modes = min(16, Nx // 3 - 1)
    x = np.linspace(0.0, 1.0, Nx, endpoint=False)
    n1 = np.zeros(Nx)
    n2 = np.zeros(Nx)
    norm_fac = math.sqrt(max(1, n_noise_modes))
    for m in range(1, n_noise_modes + 1):
        a1m, b1m = rng.standard_normal(2)
        a2m, b2m = rng.standard_normal(2)
        ang = 2.0 * np.pi * m * x
        n1 += (a1m * np.cos(ang) + b1m * np.sin(ang)) / norm_fac
        n2 += (a2m * np.cos(ang) + b2m * np.sin(ang)) / norm_fac

    if init_type == "steady_n1":
        base1, base2 = _solve_steady_n1_profile(a, a_Delta, kappa, rho, alpha, ell, p1_bar, p2_bar, Nx)
        p1_0 = base1 + noise_amp * n1
        p2_0 = base2 + noise_amp * n2
    else:
        p1_0 = p1_bar + noise_amp * n1
        p2_0 = p2_bar + noise_amp * n2

    p1_0 += p1_bar - float(np.mean(p1_0))
    p2_0 += p2_bar - float(np.mean(p2_0))

    q = 2.0 * np.pi * np.fft.rfftfreq(Nx, d=1.0 / Nx) / ell
    q2 = q * q
    dealias = np.ones_like(q)
    dealias[int(Nx // 3) :] = 0.0

    L11 = -q2 * (q2 + a + S)
    L12 = q2 * (rho + alpha)
    L21 = Q * q2 * (rho - alpha)
    L22 = -Q * q2 * (kappa * q2 + a + a_Delta + S)

    def invert_2x2(c):
        A11 = 1.0 - c * L11
        A12 = -c * L12
        A21 = -c * L21
        A22 = 1.0 - c * L22
        det = A11 * A22 - A12 * A21
        return A22 / det, -A12 / det, -A21 / det, A11 / det

    I1_11, I1_12, I1_21, I1_22 = invert_2x2(dt)
    I2_11, I2_12, I2_21, I2_22 = invert_2x2(0.5 * dt)

    p1h = np.fft.rfft(p1_0)
    p2h = np.fft.rfft(p2_0)
    p1h[0] = p1_bar * Nx
    p2h[0] = p2_bar * Nx

    def nonlin(h1, h2):
        p1 = np.fft.irfft(h1, Nx)
        p2 = np.fft.irfft(h2, Nx)
        return (
            -q2 * np.fft.rfft(p1**3 - S * p1) * dealias,
            -Q * q2 * np.fft.rfft(p2**3 - S * p2) * dealias,
        )

    nsteps = max(2, int(round(t_end / dt)))
    n1_prev, n2_prev = nonlin(p1h, p2h)
    r1 = p1h + dt * n1_prev
    r2 = p2h + dt * n2_prev
    p1h = I1_11 * r1 + I1_12 * r2
    p2h = I1_21 * r1 + I1_22 * r2

    record_start = max(1, int(0.8 * nsteps))
    record_stride = max(1, (nsteps - record_start) // 100)
    phase_hist = []
    norm_hist = []

    for step in range(1, nsteps):
        nl1, nl2 = nonlin(p1h, p2h)
        r1 = p1h + 0.5 * dt * (L11 * p1h + L12 * p2h) + dt * (1.5 * nl1 - 0.5 * n1_prev)
        r2 = p2h + 0.5 * dt * (L21 * p1h + L22 * p2h) + dt * (1.5 * nl2 - 0.5 * n2_prev)
        p1h = I2_11 * r1 + I2_12 * r2
        p2h = I2_21 * r1 + I2_22 * r2
        n1_prev, n2_prev = nl1, nl2
        if step >= record_start and (step - record_start) % record_stride == 0:
            phase_hist.append((step * dt, float(np.angle(p1h[1]))))
            p1t = np.fft.irfft(p1h, Nx)
            p2t = np.fft.irfft(p2h, Nx)
            norm_hist.append(math.sqrt(float(np.mean((p1t - p1_bar) ** 2 + (p2t - p2_bar) ** 2))))

    p1 = np.fft.irfft(p1h, Nx)
    p2 = np.fft.irfft(p2h, Nx)
    amps = np.abs(p1h[1 : min(12, Nx // 3)])
    dom_n = int(np.argmax(amps) + 1)
    norm_end = math.sqrt(float(np.mean((p1 - p1_bar) ** 2 + (p2 - p2_bar) ** 2)))
    norm_std = float(np.std(norm_hist)) if norm_hist else 0.0
    if len(phase_hist) > 1:
        ts = np.array([z[0] for z in phase_hist])
        phs = np.unwrap([z[1] for z in phase_hist])
        v_drift = abs(float((phs[-1] - phs[0]) / (2.0 * np.pi * (ts[-1] - ts[0]))))
    else:
        v_drift = 0.0

    return {
        "n_dominant_final": float(dom_n),
        "norm_final": norm_end,
        "norm_osc_std": norm_std,
        "drift_speed_final": v_drift,
        "phi1_min_final": float(np.min(p1)),
        "phi1_max_final": float(np.max(p1)),
        "phi1_mean_final": float(np.mean(p1)),
        "phi2_mean_final": float(np.mean(p2)),
    }
