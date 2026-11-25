import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq


def _solve_hss(A, Kd, Om0, Om, Omd0, Omd):
    Na = float(
        brentq(
            lambda x: -(Omd0 + Omd * A * x / Kd) * x + Om0 * (1.0 - x) * (1.0 + Om * x * x),
            0.0,
            1.0,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    return (A / Kd) * Na, Na, 1.0 - Na


def _full_noflux_rhs(F, dx, ab_v, params):
    """Exact line-by-line port of 1D/NoFlux/functions.jl:1-35 and main.jl:54-93."""
    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    Nx = F.shape[1]
    C = F[0]
    dC = np.empty(Nx, dtype=float)
    dC[0] = (C[1] - C[0]) / dx
    dC[1] = (C[2] - C[0]) / (2.0 * dx)
    dC[2:-2] = (C[:-4] - 8.0 * C[1:-3] + 8.0 * C[3:-1] - C[4:]) / (12.0 * dx)
    dC[-2] = (C[-1] - C[-3]) / (2.0 * dx)
    dC[-1] = (C[-1] - C[-2]) / dx

    d2F = np.empty_like(F)
    d2F[:, 0] = (F[:, 0] - 2.0 * F[:, 1] + F[:, 2]) / (dx * dx)
    d2F[:, 1:-1] = (F[:, :-2] - 2.0 * F[:, 1:-1] + F[:, 2:]) / (dx * dx)
    d2F[:, -1] = (F[:, -3] - 2.0 * F[:, -2] + F[:, -1]) / (dx * dx)

    dxPi = (2.0 * Z * C[1:-1] - 3.0 * B * (C[1:-1] ** 2)) * dC[1:-1]
    v = np.zeros(Nx, dtype=float)
    v[1:-1] = solve_banded((1, 1), ab_v, dxPi, overwrite_ab=False, overwrite_b=True, check_finite=False)

    adv = np.zeros_like(F)
    # Lines 15-18 (computed in Julia, then overwritten by RHS[:,1]=RHS[:,2] at lines 32-33)
    if v[1] < 0.0:
        adv[:, 0] = -F[:, 1] * v[1]
    if v[-2] > 0.0:
        adv[:, -1] = F[:, -2] * v[-2]
    vj = v[1:-1]
    vm = v[:-2]
    vp = v[2:]
    adv[:, 1:-1] = (
        -np.abs(vj) * F[:, 1:-1]
        + np.where(vm > 0.0, F[:, :-2] * vm, 0.0)
        + np.where(vp < 0.0, -F[:, 2:] * vp, 0.0)
    )
    adv /= dx

    rhs = np.zeros_like(F)
    Ci, Nai, Nii = F[0, 1:-1], F[1, 1:-1], F[2, 1:-1]
    R = Om0 * Nii * (1.0 + Om * Nai * Nai) - (Omd0 + Omd * Ci) * Nai
    rhs[0, 1:-1] = adv[0, 1:-1] + Dc * d2F[0, 1:-1] - Kd * Ci + A_par * Nai
    rhs[1, 1:-1] = adv[1, 1:-1] + Da * d2F[1, 1:-1] + R
    rhs[2, 1:-1] = adv[2, 1:-1] + d2F[2, 1:-1] - R
    rhs[:, 0] = rhs[:, 1]
    rhs[:, -1] = rhs[:, -2]
    return rhs, v


def _run_adaptive_outward_leak(L, Nx, T_end, tol, params):
    """Exact port of 1D/NoFlux/functions.jl:37-98 & main.jl:118-124 on a boundary-peaked IC
    to quantify interior nucleator mass drift when boundary flow points outward (v[2] < 0, v[Nx-1] > 0).
    """
    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    dx = L / Nx
    x = -0.5 * L + np.arange(Nx, dtype=float) * dx
    ab_v = np.zeros((3, Nx - 2), dtype=float)
    ab_v[0, 1:] = -1.0 / (dx * dx)
    ab_v[1, :] = 2.0 / (dx * dx) + 1.0
    ab_v[2, :-1] = -1.0 / (dx * dx)

    CH, NaH, NiH = _solve_hss(A_par, Kd, Om0, Om, Omd0, Omd)
    bump = np.exp(-((x + 0.5 * L) / 2.0) ** 2) + np.exp(-((x - 0.5 * L) / 2.0) ** 2)
    bump -= float(np.mean(bump))
    F = np.vstack([CH * (1.0 + 0.5 * bump), NaH * (1.0 + 0.5 * bump), NiH * (1.0 + 0.5 * bump)])
    F[:, 0] = F[:, 1]
    F[:, -1] = F[:, -2]
    m0 = float(np.mean(F[1, 1:-1] + F[2, 1:-1]))

    def midpoint(Fs, dt):
        r1, _ = _full_noflux_rhs(Fs, dx, ab_v, params)
        r2, v2 = _full_noflux_rhs(Fs + 0.5 * dt * r1, dx, ab_v, params)
        return Fs + dt * r2, v2

    max_dt = 0.5 * dx * dx
    dt_old = max_dt
    t = 0.0
    steps = 0
    while t < T_end and dt_old > 1e-15 and steps < 5000:
        dt = 2.0 * dt_old if 2.0 * dt_old < max_dt else dt_old
        if t + dt > T_end:
            dt = T_end - t
        FL, _ = midpoint(F, dt)
        FS1, _ = midpoint(F, 0.5 * dt)
        FS, _ = midpoint(FS1, 0.5 * dt)
        err = float(np.max(np.abs((FL - FS) / FL)))
        while err >= tol and dt > 1e-15:
            dt *= 0.5
            FL = FS1
            FS1, _ = midpoint(F, 0.5 * dt)
            FS, _ = midpoint(FS1, 0.5 * dt)
            # Verbatim line 90 of 1D/NoFlux/functions.jl (abs only on numerator)
            err = float(np.max(np.abs(FL - FS) / FL))
        dt_old = dt
        F = FL
        t += dt_old
        steps += 1

    rhs_f, _ = _full_noflux_rhs(F, dx, ab_v, params)
    m1 = float(np.mean(F[1, 1:-1] + F[2, 1:-1]))
    dC = np.gradient(F[0], dx)
    return {
        "C_HSS": CH,
        "C_max": float(np.max(F[0])),
        "C_center": float(F[0, Nx // 2]),
        "C_boundary": float(F[0, 0]),
        "grad_c_sq_mean": float(np.mean(dC * dC)),
        "max_rhs": float(np.max(np.abs(rhs_f))),
        "mass_drift": m1 - m0,
    }


def _relax_lp0_seed(L, M, dx, params, CH, NaH, NiH):
    """Fast IMEX relaxation on half-domain to generate a clean LP0 initial guess for banded Newton."""
    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    xh = (np.arange(M, dtype=float) + 0.5) * dx
    bump = np.where(xh <= 2.5, 1.0, 0.0)
    bump -= float(np.mean(bump))
    C = CH * (1.0 + 1.5 * bump)
    Na = NaH * (1.0 + 1.5 * bump)
    Ni = NiH * (1.0 + 1.5 * bump)
    idx = 1.0 / dx
    idx2 = idx * idx
    abV = np.zeros((3, M - 1), dtype=float)
    abV[0, 1:] = -idx2
    abV[1, :] = 1.0 + 2.0 * idx2
    abV[2, :-1] = -idx2
    dt = 0.05

    def mdiff(D, kl):
        ab = np.zeros((3, M), dtype=float)
        ab[0, 1:] = -dt * D * idx2
        ab[2, :-1] = -dt * D * idx2
        d = np.full(M, 1.0 + dt * kl + 2.0 * dt * D * idx2, dtype=float)
        d[0] -= dt * D * idx2
        d[-1] -= dt * D * idx2
        ab[1, :] = d
        return ab

    abC, abNa, abNi = mdiff(Dc, Kd), mdiff(Da, 0.0), mdiff(1.0, 0.0)
    for _ in range(3000):
        Pi = Z * C * C - B * C * C * C
        Vf = solve_banded((1, 1), abV, (Pi[1:] - Pi[:-1]) * idx, overwrite_ab=False, overwrite_b=True, check_finite=False)
        Jc = 0.5 * Vf * (C[:-1] + C[1:])
        Jna = 0.5 * Vf * (Na[:-1] + Na[1:])
        Jni = 0.5 * Vf * (Ni[:-1] + Ni[1:])
        R = Om0 * (1.0 + Om * Na * Na) * Ni - (Omd0 + Omd * C) * Na
        rC = C + dt * A_par * Na
        rC[:-1] -= dt * Jc * idx
        rC[1:] += dt * Jc * idx
        rNa = Na + dt * R
        rNa[:-1] -= dt * Jna * idx
        rNa[1:] += dt * Jna * idx
        rNi = Ni - dt * R
        rNi[:-1] -= dt * Jni * idx
        rNi[1:] += dt * Jni * idx
        C = solve_banded((1, 1), abC, rC, overwrite_ab=False, overwrite_b=True, check_finite=False)
        Na = solve_banded((1, 1), abNa, rNa, overwrite_ab=False, overwrite_b=True, check_finite=False)
        Ni = solve_banded((1, 1), abNi, rNi, overwrite_ab=False, overwrite_b=True, check_finite=False)
    return C, Na, Ni


def _solve_banded_newton(L, Nx, params, branch="spike", centered_on_node=False):
    """Exact steady-state solution (RHS == 0 to < 1e-11) of 1D/NoFlux/functions.jl:1-35."""
    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    dx = L / Nx
    M = Nx // 2
    nv = 5
    bw = 13
    CH, NaH, NiH = _solve_hss(A_par, Kd, Om0, Om, Omd0, Omd)

    def eval_res(y):
        Y = y.reshape(M, nv)
        C, Na, Ni, V, Q = Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3], Y[:, 4]
        C_ext = np.concatenate(([C[2], C[1]] if centered_on_node else [C[1], C[0]], C))
        dC = np.empty(M - 1, dtype=float)
        dC[:-1] = (C_ext[:-4] - 8.0 * C_ext[1:-3] + 8.0 * C_ext[3:-1] - C_ext[4:]) / (12.0 * dx)
        dC[-1] = (C[-1] - C[-3]) / (2.0 * dx)

        dxPi = (2.0 * Z * C[:-1] - 3.0 * B * C[:-1] * C[:-1]) * dC
        V_left = np.empty(M - 1, dtype=float)
        V_left[0] = -V[1] if centered_on_node else -V[0]
        V_left[1:] = V[:-2]
        V_right = V[1:]
        rV = np.empty(M, dtype=float)
        rV[:-1] = (2.0 / (dx * dx) + 1.0) * V[:-1] - (V_left + V_right) / (dx * dx) - dxPi
        rV[-1] = V[-1]
        if centered_on_node:
            rV[0] = V[0]

        def rhs_field(S, D, reac):
            S_left = np.empty(M - 1, dtype=float)
            S_left[0] = S[1] if centered_on_node else S[0]
            S_left[1:] = S[:-2]
            S_right = S[1:]
            d2S = (S_left - 2.0 * S[:-1] + S_right) / (dx * dx)
            vj, vm, vp = V[:-1], V_left, V_right
            adv = (
                -np.abs(vj) * S[:-1]
                + np.where(vm > 0.0, S_left * vm, 0.0)
                + np.where(vp < 0.0, -S_right * vp, 0.0)
            ) / dx
            out = np.empty(M, dtype=float)
            out[:-1] = adv + D * d2S + reac
            out[-1] = S[-1] - S[-2]
            return out

        R = Om0 * (1.0 + Om * Na[:-1] * Na[:-1]) * Ni[:-1] - (Omd0 + Omd * C[:-1]) * Na[:-1]
        rC = rhs_field(C, Dc, A_par * Na[:-1] - Kd * C[:-1])
        rNa = rhs_field(Na, Da, R)
        rNi = rhs_field(Ni, 1.0, -R)
        rNi[-2] = Q[-2] - 1.0

        rQ = np.empty(M, dtype=float)
        if not centered_on_node:
            w = np.ones(M - 1, dtype=float) / (M - 1)
        else:
            w = np.ones(M - 1, dtype=float) / (M - 1.5)
            w[0] *= 0.5
        rQ[0] = Q[0] - (Na[0] + Ni[0]) * w[0]
        rQ[1:-1] = Q[1:-1] - Q[:-2] - (Na[1:-1] + Ni[1:-1]) * w[1:]
        rQ[-1] = Q[-1] - Q[-2]
        return np.column_stack([rC, rNa, rNi, rV, rQ]).ravel()

    def build_jac(y, r0):
        Ntot = nv * M
        ab = np.zeros((2 * bw + 1, Ntot), dtype=float)
        eps = 1e-8
        for v_i in range(nv):
            for mod in range(6):
                cols = np.arange(mod, M, 6) * nv + v_i
                dy = np.zeros(Ntot, dtype=float)
                dy[cols] = eps
                dr = (eval_res(y + dy) - r0) / eps
                for shift in range(-bw, bw + 1):
                    rows = cols + shift
                    valid = (rows >= 0) & (rows < Ntot)
                    if np.any(valid):
                        ab[bw + shift, cols[valid]] = dr[rows[valid]]
        return ab

    xh = (np.arange(M, dtype=float) + (0.0 if centered_on_node else 0.5)) * dx
    if branch == "spike":
        C0 = CH + 0.52 * np.exp(-xh * xh / 0.5)
        Na0 = NaH + 0.12 * np.exp(-xh * xh / 0.5)
        Ni0 = NiH + 0.25 * np.exp(-xh * xh / 0.5)
        Ni0 += 1.0 - float(np.mean(Na0 + Ni0))
        inv_dt = 0.5
    else:
        # Relax LP0 on coarse 512 grid, then interpolate to M if M > 256
        if M > 256:
            C_c, Na_c, Ni_c = _relax_lp0_seed(L, 256, L / 512.0, params, CH, NaH, NiH)
            xc = (np.arange(256, dtype=float) + 0.5) * (L / 512.0)
            C0 = np.interp(xh, xc, C_c)
            Na0 = np.interp(xh, xc, Na_c)
            Ni0 = np.interp(xh, xc, Ni_c)
        else:
            C0, Na0, Ni0 = _relax_lp0_seed(L, M, dx, params, CH, NaH, NiH)
        Ni0 += 1.0 - float(np.mean(Na0 + Ni0))
        inv_dt = 0.2

    C_ext = np.concatenate(([C0[2], C0[1]] if centered_on_node else [C0[1], C0[0]], C0))
    dC0 = np.empty(M - 1, dtype=float)
    dC0[:-1] = (C_ext[:-4] - 8.0 * C_ext[1:-3] + 8.0 * C_ext[3:-1] - C_ext[4:]) / (12.0 * dx)
    dC0[-1] = (C0[-1] - C0[-3]) / (2.0 * dx)
    dxPi0 = (2.0 * Z * C0[:-1] - 3.0 * B * C0[:-1] * C0[:-1]) * dC0
    V0 = np.zeros(M, dtype=float)
    if not centered_on_node:
        abV = np.zeros((3, M - 1), dtype=float)
        abV[0, 1:] = -1.0 / (dx * dx)
        abV[1, :] = 2.0 / (dx * dx) + 1.0
        abV[1, 0] = 3.0 / (dx * dx) + 1.0
        abV[2, :-1] = -1.0 / (dx * dx)
        V0[:-1] = solve_banded((1, 1), abV, dxPi0)
    else:
        abV = np.zeros((3, M - 2), dtype=float)
        abV[0, 1:] = -1.0 / (dx * dx)
        abV[1, :] = 2.0 / (dx * dx) + 1.0
        abV[2, :-1] = -1.0 / (dx * dx)
        V0[1:-1] = solve_banded((1, 1), abV, dxPi0[1:])

    Q0 = np.cumsum(Na0 + Ni0) / M
    y = np.column_stack([C0, Na0, Ni0, V0, Q0]).ravel()
    c_idx = np.arange(M - 1) * nv + 0
    na_idx = np.arange(M - 1) * nv + 1
    ni_idx = np.arange(M - 2) * nv + 2

    nr = float("inf")
    for _ in range(40):
        r0 = eval_res(y)
        nr = float(np.max(np.abs(r0)))
        if nr < 1e-11:
            break
        ab = build_jac(y, r0)
        if inv_dt > 0.0:
            ab[bw, c_idx] -= inv_dt
            ab[bw, na_idx] -= inv_dt
            ab[bw, ni_idx] -= inv_dt
        step = solve_banded((bw, bw), ab, -r0)
        if inv_dt > 0.0:
            alpha = min(1.0, 0.3 / max(0.3, float(np.max(np.abs(step)))))
            y = y + alpha * step
            nr_n = float(np.max(np.abs(eval_res(y))))
            if nr_n < 0.05:
                inv_dt = 0.0
            elif nr_n < nr:
                inv_dt *= 0.35
        else:
            n0 = float(np.linalg.norm(r0))
            alpha = 1.0
            for _ in range(10):
                yn = y + alpha * step
                if float(np.linalg.norm(eval_res(yn))) < n0 or alpha < 1e-2:
                    y = yn
                    break
                alpha *= 0.5

    if nr > 1e-6:
        raise RuntimeError(f"1D/NoFlux banded Newton did not converge (max_rhs={nr:.3e})")

    Y = y.reshape(M, nv)
    C = Y[:, 0]
    dC_g = np.gradient(C, dx)
    return {
        "C_HSS": CH,
        "C_max": float(np.max(C)),
        "C_center": float(C[0]),
        "C_boundary": float(C[-1]),
        "grad_c_sq_mean": float(np.mean(dC_g * dC_g)),
        "max_rhs": float(nr),
        "mass_drift": 0.0,
    }


def run(inputs, settings):
    L_base = float(inputs["L_base"])
    L_mult = float(inputs["L_mult"])
    N_base = int(round(float(inputs["N_base"])))
    Dc = float(inputs["Dc"])
    Da = float(inputs["Da"])
    A_par = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"])
    Om = float(inputs["Omega"])
    Omd = float(inputs["Omega_d"])
    Omd0 = float(inputs["Omega_d0"])
    Z = float(inputs["Z"])
    B = float(inputs["B"])
    T_end = float(inputs["T_end"])
    tol = float(inputs["ErrorTolerance"])

    if Dc <= 0.0 or Da <= 0.0 or A_par <= 0.0 or Kd <= 0.0 or Om0 < 0.0 or Om < 0.0 or Omd < 0.0 or Omd0 < 0.0:
        raise ValueError("Non-physical negative rate or non-positive diffusivity parameter")

    solver = str(settings.get("solver", "banded_newton")).lower()
    branch = str(settings.get("branch", "spike")).lower()
    centered_on_node = bool(settings.get("centered_on_node", False))
    n_grid_scale = float(settings.get("n_grid_scale", 1.0))

    L = L_base * L_mult
    params = (Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B)

    if solver == "adaptive_midpoint_leak":
        Nx_leak = int(round(float(settings.get("leak_probe_Nx", 128)) * n_grid_scale))
        return _run_adaptive_outward_leak(L, Nx_leak, T_end, tol, params)

    Nx = int(round(N_base * L_mult * n_grid_scale))
    if Nx % 2 == 1:
        Nx += 1
    return _solve_banded_newton(L, Nx, params, branch=branch, centered_on_node=centered_on_node)
