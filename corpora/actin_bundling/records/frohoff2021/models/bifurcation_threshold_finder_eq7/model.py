import math
import numpy as np
from scipy.optimize import brentq, root


def _warm_start(a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, amp1, amp2):
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


def _make_F(a_Delta, kappa, rho, alpha, ell, p1_bar, p2_bar, n_per, M):
    N = 4 * M
    theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    ms = np.arange(1, M + 1)
    C = np.cos(np.outer(theta, ms))
    k_sq = (2.0 * np.pi * n_per * ms / ell) ** 2

    def F(u, a):
        c1, c2 = u[:M], u[M:]
        p1 = p1_bar + C @ c1
        p2 = p2_bar + C @ c2
        f1p = a * p1 + p1**3 - (rho + alpha) * p2
        f2p = (a + a_Delta) * p2 + p2**3 - (rho - alpha) * p1
        r1 = k_sq * c1 + (2.0 / N) * (C.T @ f1p)
        r2 = kappa * k_sq * c2 + (2.0 / N) * (C.T @ f2p)
        return np.concatenate([r1, r2])

    return F


def _bloch_evals(c1, c2, a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, N_cell):
    M = len(c1)
    y = np.linspace(0.0, 1.0, N_cell, endpoint=False)
    ms = np.arange(1, M + 1)
    C = np.cos(np.outer(2.0 * np.pi * y, ms))
    p1 = p1_bar + C @ c1
    p2 = p2_bar + C @ c2
    F_mat = np.fft.fft(np.eye(N_cell), axis=0)
    Finv = np.fft.ifft(np.eye(N_cell), axis=0)
    m_idx = np.fft.fftfreq(N_cell, d=1.0 / N_cell)
    all_ev = []
    for j in range(n_per):
        q = 2.0 * np.pi * (n_per * m_idx + j) / ell
        D2 = Finv @ np.diag(-(q * q)) @ F_mat
        L11 = -D2 + np.diag(a + 3.0 * p1 * p1)
        L12 = -(rho + alpha) * np.eye(N_cell)
        L21 = -(rho - alpha) * np.eye(N_cell)
        L22 = -kappa * D2 + np.diag(a + a_Delta + 3.0 * p2 * p2)
        J = np.block([[D2 @ L11, D2 @ L12], [Q * D2 @ L21, Q * D2 @ L22]])
        ev = np.linalg.eigvals(J)
        if j == 0:
            order = np.argsort(np.abs(ev))
            ev = ev[order[3:]]
        all_ev.append(ev)
    evals = np.concatenate(all_ev)
    return np.sort(evals.real)[::-1]


def run(inputs, settings):
    a_start = inputs["a_start"]
    a_end = inputs["a_end"]
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

    mode = str(settings.get("mode", "stability_1st"))
    M = int(settings.get("n_modes", 20))
    N_cell = max(int(settings.get("n_cell", 48)), 2 * M + 4)
    ds_max = float(settings.get("ds_max", 0.03))
    tol = float(settings.get("tol", 1e-11))
    warm_start = bool(settings.get("warm_start", False))

    F = _make_F(a_Delta, kappa, rho, alpha, ell, p1_bar, p2_bar, n_per, M)

    def get_u0(a_val):
        if warm_start:
            return _warm_start(a_val, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, amp1, amp2)
        u0 = np.zeros(2 * M)
        u0[0] = amp1
        u0[M] = amp2
        return u0

    if mode in ("stability_1st", "stability_2nd", "drift_pitchfork"):
        u_init = root(lambda u: F(u, a_start), get_u0(a_start), method="hybr", tol=tol).x
        cache_u = [u_init]

        def solve_at(a_val):
            sol = root(lambda u: F(u, a_val), cache_u[0], method="hybr", tol=tol)
            if not sol.success or math.sqrt(0.5 * float(np.sum(sol.x**2))) < 1e-4:
                sol = root(lambda u: F(u, a_val), get_u0(a_val), method="hybr", tol=tol)
            cache_u[0] = sol.x
            return sol.x

        def obj(a_val):
            u = solve_at(a_val)
            c1, c2 = u[:M], u[M:]
            if mode == "drift_pitchfork":
                return 0.5 * float(np.sum(c1 * c1) + (1.0 / Q) * ((rho + alpha) / (rho - alpha)) * np.sum(c2 * c2))
            ev = _bloch_evals(c1, c2, a_val, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, N_cell)
            return float(ev[0] if mode == "stability_1st" else ev[1])

        a_lo, a_hi = min(a_start, a_end), max(a_start, a_end)
        a_crit = float(brentq(obj, a_lo, a_hi, xtol=1e-9, rtol=1e-9))
        u_c = solve_at(a_crit)
        norm_c = math.sqrt(0.5 * float(np.sum(u_c**2)))
        return {"a_crit": a_crit, "a_crit_2": None, "norm_at_crit": norm_c}

    if mode == "fold":
        u0 = get_u0(a_start)
        s0 = root(lambda u: F(u, a_start), u0, method="hybr", tol=tol)
        u_prev = s0.x
        step_dir = 1.0 if a_end >= a_start else -1.0
        da0 = step_dir * min(abs(a_end - a_start) * 0.05, 0.5 * ds_max)
        s1 = root(lambda u: F(u, a_start + da0), u_prev, method="hybr", tol=tol)
        z_prev = np.append(u_prev, a_start)
        z_curr = np.append(s1.x, a_start + da0)
        folds = []
        a_lo_lim = min(a_start, a_end) - 0.5
        a_hi_lim = max(a_start, a_end) + 0.5
        for _ in range(800):
            t = z_curr - z_prev
            nrm = float(np.linalg.norm(t))
            if nrm < 1e-13:
                break
            t = t / nrm
            ds = min(ds_max, max(0.002, nrm))
            z_pred = z_curr + ds * t

            def aug(z):
                return np.append(F(z[:-1], z[-1]), np.dot(z - z_curr, t) - ds)

            sol = root(aug, z_pred, method="hybr", tol=tol)
            if not sol.success:
                ds *= 0.5
                z_pred = z_curr + ds * t
                sol = root(aug, z_pred, method="hybr", tol=tol)
                if not sol.success:
                    break
            z_next = sol.x
            da_prev = z_curr[-1] - z_prev[-1]
            da_next = z_next[-1] - z_curr[-1]
            if da_prev * da_next < 0.0:
                # Quadratic interpolation on arclength (0, s1, s2)
                s1_len = float(np.linalg.norm(z_curr - z_prev))
                s2_len = float(np.linalg.norm(z_next - z_curr))
                a0, a1, a2 = float(z_prev[-1]), float(z_curr[-1]), float(z_next[-1])
                # Fit quadratic a(s) through (-s1, a0), (0, a1), (s2, a2)
                denom = s1_len * s2_len * (s1_len + s2_len)
                A_q = (s1_len * (a2 - a1) + s2_len * (a0 - a1)) / denom
                B_q = (s1_len**2 * (a2 - a1) - s2_len**2 * (a0 - a1)) / denom
                s_star = -B_q / (2.0 * A_q) if abs(A_q) > 1e-15 else 0.0
                a_fold = a1 + B_q * s_star + A_q * s_star**2
                norm_fold = math.sqrt(0.5 * float(np.sum(z_curr[:-1] ** 2)))
                folds.append((a_fold, norm_fold))
            z_prev, z_curr = z_curr, z_next
            norm_now = math.sqrt(0.5 * float(np.sum(z_curr[:-1] ** 2)))
            if norm_now < 0.02 or z_curr[-1] < a_lo_lim or z_curr[-1] > a_hi_lim:
                break
            if len(folds) >= 2:
                break
        if not folds:
            return {"a_crit": None, "a_crit_2": None, "norm_at_crit": None}
        return {
            "a_crit": float(folds[0][0]),
            "a_crit_2": float(folds[1][0]) if len(folds) >= 2 else None,
            "norm_at_crit": float(folds[0][1]),
        }

    raise ValueError(f"Unknown mode: {mode}")
