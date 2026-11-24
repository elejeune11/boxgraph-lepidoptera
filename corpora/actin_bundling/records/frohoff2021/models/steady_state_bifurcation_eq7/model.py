import math
import numpy as np
from scipy.optimize import root


def _warm_start_profile(a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, amp1, amp2):
    N = 4 * M
    theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    p1 = p1_bar + amp1 * np.cos(theta)
    p2 = p2_bar + amp2 * np.cos(theta)
    q = 2.0 * np.pi * n_per * np.fft.rfftfreq(N, d=1.0 / N) / ell
    q2 = q * q
    S = 6.0
    dt = 0.05
    L11 = -q2 * (q2 + a + S)
    L12 = q2 * (rho + alpha)
    L21 = Q * q2 * (rho - alpha)
    L22 = -Q * q2 * (kappa * q2 + a + a_Delta + S)
    A11 = 1.0 - dt * L11
    A12 = -dt * L12
    A21 = -dt * L21
    A22 = 1.0 - dt * L22
    det = A11 * A22 - A12 * A21
    I11, I12, I21, I22 = A22 / det, -A12 / det, -A21 / det, A11 / det
    h1 = np.fft.rfft(p1)
    h2 = np.fft.rfft(p2)
    h1[0] = p1_bar * N
    h2[0] = p2_bar * N
    for _ in range(400):
        p1 = np.fft.irfft(h1, N)
        p2 = np.fft.irfft(h2, N)
        n1 = -q2 * np.fft.rfft(p1**3 - S * p1)
        n2 = -Q * q2 * np.fft.rfft(p2**3 - S * p2)
        r1 = h1 + dt * n1
        r2 = h2 + dt * n2
        h1 = I11 * r1 + I12 * r2
        h2 = I21 * r1 + I22 * r2
    p1 = np.fft.irfft(h1, N)
    p2 = np.fft.irfft(h2, N)
    idx = int(np.argmax(p1))
    p1 = np.roll(p1, -idx)
    p2 = np.roll(p2, -idx)
    c1 = (2.0 / N) * np.real(np.fft.rfft(p1)[1 : M + 1])
    c2 = (2.0 / N) * np.real(np.fft.rfft(p2)[1 : M + 1])
    return np.concatenate([c1, c2])


def solve_even_steady(a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, tol, amp1, amp2, warm_start=False, init_u=None):
    N = 4 * M
    theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    ms = np.arange(1, M + 1)
    C = np.cos(np.outer(theta, ms))
    k_sq = (2.0 * np.pi * n_per * ms / ell) ** 2

    def res(u):
        c1, c2 = u[:M], u[M:]
        p1 = p1_bar + C @ c1
        p2 = p2_bar + C @ c2
        f1p = a * p1 + p1**3 - (rho + alpha) * p2
        f2p = (a + a_Delta) * p2 + p2**3 - (rho - alpha) * p1
        r1 = k_sq * c1 + (2.0 / N) * (C.T @ f1p)
        r2 = kappa * k_sq * c2 + (2.0 / N) * (C.T @ f2p)
        return np.concatenate([r1, r2])

    if init_u is not None:
        u0 = init_u
    elif warm_start:
        u0 = _warm_start_profile(a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, amp1, amp2)
    else:
        u0 = np.zeros(2 * M)
        u0[0] = amp1
        u0[M] = amp2

    sol = root(res, u0, method="hybr", tol=tol)
    return sol.x, bool(sol.success and np.max(np.abs(sol.fun)) < max(1e-7, 1e3 * tol))


def bloch_eigenvalues(c1, c2, a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, N_cell):
    M = len(c1)
    y = np.linspace(0.0, 1.0, N_cell, endpoint=False)
    ms = np.arange(1, M + 1)
    C = np.cos(np.outer(2.0 * np.pi * y, ms))
    p1 = p1_bar + C @ c1
    p2 = p2_bar + C @ c2
    F = np.fft.fft(np.eye(N_cell), axis=0)
    Finv = np.fft.ifft(np.eye(N_cell), axis=0)
    m_idx = np.fft.fftfreq(N_cell, d=1.0 / N_cell)
    all_ev = []
    for j in range(n_per):
        q = 2.0 * np.pi * (n_per * m_idx + j) / ell
        D2 = Finv @ np.diag(-(q * q)) @ F
        L11 = -D2 + np.diag(a + 3.0 * p1 * p1)
        L12 = -(rho + alpha) * np.eye(N_cell)
        L21 = -(rho - alpha) * np.eye(N_cell)
        L22 = -kappa * D2 + np.diag(a + a_Delta + 3.0 * p2 * p2)
        J = np.block([[D2 @ L11, D2 @ L12], [Q * D2 @ L21, Q * D2 @ L22]])
        ev = np.linalg.eigvals(J)
        if j == 0:
            order = np.argsort(np.abs(ev))
            ev = ev[order[3:]]  # drop 2 mass-conservation + 1 translation zero modes
        all_ev.append(ev)
    evals = np.concatenate(all_ev)
    idx = np.argsort(evals.real)[::-1]
    return evals[idx], p1, p2


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
    n_per = max(1, int(round(inputs["n_period"])))
    amp1 = inputs["init_amp_1"]
    amp2 = inputs["init_amp_2"]

    M = int(settings.get("n_modes", 20))
    N_cell = max(int(settings.get("n_cell", 48)), 2 * M + 4)
    tol = float(settings.get("tol", 1e-12))
    warm_start = bool(settings.get("warm_start", False))

    u, ok = solve_even_steady(
        a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, tol, amp1, amp2, warm_start=warm_start
    )
    if not ok:
        return {k: None for k in ("norm_dphi", "lam_max", "n_unstable", "dp_integral", "phi1_min", "phi1_max", "phi2_min", "phi2_max")}

    c1, c2 = u[:M], u[M:]
    evals, p1, p2 = bloch_eigenvalues(c1, c2, a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, N_cell)
    norm_dphi = math.sqrt(0.5 * float(np.sum(c1 * c1 + c2 * c2)))
    n_unstable = int(np.sum(evals.real > 1e-8))
    if abs(rho - alpha) > 1e-15:
        dp_int = 0.5 * float(np.sum(c1 * c1) + (1.0 / Q) * ((rho + alpha) / (rho - alpha)) * np.sum(c2 * c2))
    else:
        dp_int = None

    return {
        "norm_dphi": norm_dphi,
        "lam_max": float(evals[0].real),
        "n_unstable": float(n_unstable),
        "dp_integral": dp_int,
        "phi1_min": float(np.min(p1)),
        "phi1_max": float(np.max(p1)),
        "phi2_min": float(np.min(p2)),
        "phi2_max": float(np.max(p2)),
    }
