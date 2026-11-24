import importlib.util
import math
from pathlib import Path
import numpy as np

_SS_PATH = Path(__file__).resolve().parent.parent / "steady_state_bifurcation_eq7" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_bg_steady_state_bifurcation_eq7", _SS_PATH)
_SS_MOD = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_SS_MOD)
solve_even_steady = _SS_MOD.solve_even_steady


def _relax_warm_start(a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, amp1, amp2, n_steps, dt):
    N = 4 * M
    theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    p1 = p1_bar + amp1 * np.cos(theta)
    p2 = p2_bar + amp2 * np.cos(theta)
    q = 2.0 * np.pi * n_per * np.fft.rfftfreq(N, d=1.0 / N) / ell
    q2 = q * q
    S = max(6.0, 3.0 * (abs(a) + abs(a_Delta) + 2.0))
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
    for _ in range(n_steps):
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
    # Center the even profile at the maximum of p1 (or p2 if p1 is flat)
    idx = int(np.argmax(p1)) if (np.max(p1) - np.min(p1)) > 1e-4 else int(np.argmax(p2))
    p1 = np.roll(p1, -idx)
    p2 = np.roll(p2, -idx)
    c1 = (2.0 / N) * np.real(np.fft.rfft(p1)[1 : M + 1])
    c2 = (2.0 / N) * np.real(np.fft.rfft(p2)[1 : M + 1])
    return np.concatenate([c1, c2])


def _linear_eigenvector_seed(a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, scale, branch):
    q = 2.0 * np.pi * n_per / ell
    q2 = q * q
    f1pp = a + 3.0 * p1_bar**2
    f2pp = a + a_Delta + 3.0 * p2_bar**2
    J = np.array(
        [
            [-q2 * (q2 + f1pp), q2 * (rho + alpha)],
            [Q * q2 * (rho - alpha), -Q * q2 * (kappa * q2 + f2pp)],
        ],
        dtype=float,
    )
    vals, vecs = np.linalg.eig(J)
    order = np.argsort(np.real(vals))[::-1]
    idx = order[0] if branch == "plus" else order[1]
    v = np.real(vecs[:, idx])
    v = v / max(np.linalg.norm(v), 1e-12) * scale
    u0 = np.zeros(2 * M)
    u0[0] = v[0]
    u0[M] = v[1]
    return u0


def _solve_with_jac(a, a_Delta, kappa, rho, alpha, ell, p1_bar, p2_bar, n_per, M, tol, u0):
    from scipy.optimize import root

    N = 4 * M
    theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    ms = np.arange(1, M + 1)
    C = np.cos(np.outer(theta, ms))
    k_sq = (2.0 * np.pi * n_per * ms / ell) ** 2
    I_M = np.eye(M)
    scale = 2.0 / N

    def res_and_jac(u):
        c1, c2 = u[:M], u[M:]
        p1 = p1_bar + C @ c1
        p2 = p2_bar + C @ c2
        f1p = a * p1 + p1**3 - (rho + alpha) * p2
        f2p = (a + a_Delta) * p2 + p2**3 - (rho - alpha) * p1
        r1 = k_sq * c1 + scale * (C.T @ f1p)
        r2 = kappa * k_sq * c2 + scale * (C.T @ f2p)
        d1 = a + 3.0 * p1 * p1
        d2 = a + a_Delta + 3.0 * p2 * p2
        J11 = np.diag(k_sq) + scale * (C.T * d1) @ C
        J12 = -(rho + alpha) * I_M
        J21 = -(rho - alpha) * I_M
        J22 = np.diag(kappa * k_sq) + scale * (C.T * d2) @ C
        J = np.block([[J11, J12], [J21, J22]])
        return np.concatenate([r1, r2]), J

    sol = root(res_and_jac, u0, jac=True, method="hybr", tol=tol)
    return sol.x, bool(sol.success and np.max(np.abs(sol.fun)) < max(1e-7, 1e3 * tol))


def continue_steady(
    a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, tol, amp1, amp2, settings
):
    start_a = settings.get("start_a", None)
    start_p2 = settings.get("start_phi2_bar", None)
    step_size = float(settings.get("step_size", 0.1))
    warm_steps = int(settings.get("warm_start_steps", 0))
    warm_dt = float(settings.get("warm_start_dt", 0.2))
    linear_seed = settings.get("linear_seed", None)

    a0 = a if start_a is None else float(start_a)
    p20 = p2_bar if start_p2 is None else float(start_p2)

    if warm_steps > 0:
        u0 = _relax_warm_start(a0, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p20, n_per, M, amp1, amp2, warm_steps, warm_dt)
    elif linear_seed in ("plus", "minus"):
        scale = math.hypot(amp1, amp2) * 0.25
        u0 = _linear_eigenvector_seed(a0, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p20, n_per, M, scale, linear_seed)
    else:
        u0 = np.zeros(2 * M)
        u0[0] = amp1
        u0[M] = amp2

    def _solve(ca, cp2, init_u):
        if M > 32:
            return _solve_with_jac(ca, a_Delta, kappa, rho, alpha, ell, p1_bar, cp2, n_per, M, tol, init_u)
        return solve_even_steady(ca, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, cp2, n_per, M, tol, amp1, amp2, init_u=init_u)

    u, ok = _solve(a0, p20, u0)
    if not ok or math.sqrt(0.5 * float(np.sum(u * u))) < 1e-4:
        return None, False

    def march(cur_val, target_val, is_a):
        nonlocal u
        u_prev = None
        val_prev = None
        h = step_size
        while abs(target_val - cur_val) > 1e-12:
            diff = target_val - cur_val
            step = math.copysign(min(abs(diff), h), diff)
            nxt = cur_val + step
            if u_prev is not None and abs(cur_val - val_prev) > 1e-14:
                u_pred = u + (u - u_prev) * (step / (cur_val - val_prev))
            else:
                u_pred = u
            ca = nxt if is_a else a0
            cp2 = p20 if is_a else nxt
            u_new, ok_new = _solve(ca, cp2, u_pred)
            norm_new = math.sqrt(0.5 * float(np.sum(u_new * u_new))) if ok_new else 0.0
            if ok_new and norm_new > 1e-4:
                u_prev, val_prev = u, cur_val
                u, cur_val = u_new, nxt
                h = min(step_size, h * 1.4)
            else:
                h *= 0.5
                if h < 1e-4:
                    return False
        return True

    if start_a is not None:
        if not march(a0, a, is_a=True):
            return None, False
        a0 = a
    if start_p2 is not None:
        if not march(p20, p2_bar, is_a=False):
            return None, False

    return u, True


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

    M = int(settings.get("n_modes", 24))
    tol = float(settings.get("tol", 1e-11))

    u, ok = continue_steady(
        a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, tol, amp1, amp2, settings
    )
    if not ok:
        return {k: None for k in ("norm_dphi", "phi1_min", "phi1_max", "phi2_min", "phi2_max")}

    c1, c2 = u[:M], u[M:]
    norm_dphi = math.sqrt(0.5 * float(np.sum(c1 * c1 + c2 * c2)))
    N_eval = max(64, 4 * M)
    theta = np.linspace(0.0, 2.0 * np.pi, N_eval, endpoint=False)
    ms = np.arange(1, M + 1)
    C = np.cos(np.outer(theta, ms))
    p1 = p1_bar + C @ c1
    p2 = p2_bar + C @ c2
    return {
        "norm_dphi": norm_dphi,
        "phi1_min": float(np.min(p1)),
        "phi1_max": float(np.max(p1)),
        "phi2_min": float(np.min(p2)),
        "phi2_max": float(np.max(p2)),
    }
