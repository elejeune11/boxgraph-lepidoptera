import math
import numpy as np
from scipy.fft import rfft, irfft, rfftfreq
from scipy.optimize import brentq, newton_krylov
from scipy.sparse.linalg import LinearOperator


def _solve_hss(A, Kd, Om0, Om, Omd0, Omd, nbar=1.0):
    Na = float(
        brentq(
            lambda x: -(Omd0 + Omd * A * x / Kd) * x + Om0 * (nbar - x) * (1.0 + Om * x * x),
            0.0,
            nbar,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    return (A / Kd) * Na, Na, nbar - Na


def _periodic_rhs(F, k, ik, k2, vel_op, params):
    """Exact line-by-line port of 1D/PeriodicBoundaries/functions.jl:1-18."""
    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    Nx = F.shape[1]
    d2 = irfft(rfft(F, axis=1) * (-k2), n=Nx, axis=1)
    sig = Z * (F[0] ** 2) - B * (F[0] ** 3)
    v = irfft(rfft(sig) * vel_op, n=Nx)
    dxFv = irfft(rfft(F * v, axis=1) * ik, n=Nx, axis=1)
    rhs = np.empty_like(F)
    ra = -(Omd0 + Omd * F[0]) * F[1] + Om0 * F[2] * (1.0 + Om * (F[1] ** 2))
    rhs[0] = -dxFv[0] + Dc * d2[0] - Kd * F[0] + A_par * F[1]
    rhs[1] = -dxFv[1] + Da * d2[1] + ra
    rhs[2] = -dxFv[2] + d2[2] - ra
    return rhs, v


def _run_adaptive_midpoint(L, Nx, dt_init, tol, params, adapt_steps=5):
    """Exact port of 1D/PeriodicBoundaries/{InputParameters,functions,main}.jl on ZeroMeanBump IC."""
    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    dx = L / Nx
    x = -0.5 * L + np.arange(Nx, dtype=float) * dx
    k = 2.0 * math.pi * rfftfreq(Nx, d=dx)
    ik = 1j * k
    k2 = k * k
    vel_op = ik / (1.0 + k2)

    CH, NaH, NiH = _solve_hss(A_par, Kd, Om0, Om, Omd0, Omd, 1.0)
    bump = np.where((x >= -5.0) & (x <= 5.0), 1.0, 0.0)
    # Verbatim InputParameters.jl:44 using trapezoidal integrate(x, Bump)/L
    bump = bump - float(np.trapezoid(bump, x)) / L
    F = np.vstack([CH * (1.0 + bump), NaH * (1.0 + bump), NiH * (1.0 + bump)])

    def midpoint(Fs, dt):
        r1, _ = _periodic_rhs(Fs, k, ik, k2, vel_op, params)
        r2, v2 = _periodic_rhs(Fs + 0.5 * dt * r1, k, ik, k2, vel_op, params)
        return Fs + dt * r2, v2

    dt_old = dt_init
    dt_first = dt_init
    for step_i in range(adapt_steps):
        dt = 2.0 * dt_old
        FL, _ = midpoint(F, dt)
        FS1, _ = midpoint(F, 0.5 * dt)
        FS, _ = midpoint(FS1, 0.5 * dt)
        err = float(np.max(np.abs((FL - FS) / FL)))
        while err >= tol and dt > 1e-15:
            dt *= 0.5
            FL = FS1
            FS1, _ = midpoint(F, 0.5 * dt)
            FS, _ = midpoint(FS1, 0.5 * dt)
            # Verbatim line 73 of 1D/PeriodicBoundaries/functions.jl
            err = float(np.max(np.abs(FL - FS) / FL))
        dt_old = dt
        if step_i == 0:
            dt_first = dt_old
        F = FL

    rhs_f, _ = _periodic_rhs(F, k, ik, k2, vel_op, params)
    dC = irfft(rfft(F[0]) * ik, n=Nx)
    return {
        "C_HSS": CH,
        "C_min": float(np.min(F[0])),
        "C_max": float(np.max(F[0])),
        "C_center": float(F[0, Nx // 2]),
        "grad_c_sq_mean": float(np.mean(dC * dC)),
        "mean_N": float(np.mean(F[1] + F[2])),
        "max_rhs": float(np.max(np.abs(rhs_f))),
        "dt_first_step": float(dt_first),
    }


def _run_imex_time(L, Nx, T_end, dt, params, ic_mode="bump", seed=123):
    """Time integration of exact 1D/PeriodicBoundaries/functions.jl Fourier RHS using
    semi-implicit diffusion (for travelling waves and oscillatory localized states).
    """
    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    dx = L / Nx
    x = -0.5 * L + np.arange(Nx, dtype=float) * dx
    k = 2.0 * math.pi * rfftfreq(Nx, d=dx)
    ik = 1j * k
    k2 = k * k
    vel_op = ik / (1.0 + k2)
    CH, NaH, NiH = _solve_hss(A_par, Kd, Om0, Om, Omd0, Omd, 1.0)

    if ic_mode == "noise":
        rng = np.random.default_rng(seed)
        noise = 0.01 * (-1.0 + 2.0 * rng.random(Nx))
        noise -= float(np.mean(noise))
        C = CH * (1.0 + noise)
        Na = NaH * (1.0 + noise)
        Ni = NiH * (1.0 + noise)
    else:
        bump = np.exp(-((x / 4.0) ** 2))
        bump -= float(np.mean(bump))
        C = CH * (1.0 + 1.5 * bump)
        Na = NaH * (1.0 + 1.5 * bump)
        Ni = NiH * (1.0 + 1.5 * bump)

    den_c = 1.0 + dt * (Dc * k2 + Kd)
    den_a = 1.0 + dt * Da * k2
    den_i = 1.0 + dt * k2
    nsteps = max(1, int(round(T_end / dt)))
    dt_s = T_end / nsteps
    den_c = 1.0 + dt_s * (Dc * k2 + Kd)
    den_a = 1.0 + dt_s * Da * k2
    den_i = 1.0 + dt_s * k2

    for _ in range(nsteps):
        sig = Z * C * C - B * C * C * C
        v = irfft(rfft(sig) * vel_op, n=Nx)
        ac = -irfft(rfft(C * v) * ik, n=Nx)
        aa = -irfft(rfft(Na * v) * ik, n=Nx)
        ai = -irfft(rfft(Ni * v) * ik, n=Nx)
        ra = -(Omd0 + Omd * C) * Na + Om0 * Ni * (1.0 + Om * Na * Na)
        C = irfft((rfft(C) + dt_s * rfft(ac + A_par * Na)) / den_c, n=Nx)
        Na = irfft((rfft(Na) + dt_s * rfft(aa + ra)) / den_a, n=Nx)
        Ni = irfft((rfft(Ni) + dt_s * rfft(ai - ra)) / den_i, n=Nx)

    F = np.vstack([C, Na, Ni])
    rhs_f, _ = _periodic_rhs(F, k, ik, k2, vel_op, params)
    dC = irfft(rfft(C) * ik, n=Nx)
    return {
        "C_HSS": CH,
        "C_min": float(np.min(C)),
        "C_max": float(np.max(C)),
        "C_center": float(C[Nx // 2]),
        "grad_c_sq_mean": float(np.mean(dC * dC)),
        "mean_N": float(np.mean(Na + Ni)),
        "max_rhs": float(np.max(np.abs(rhs_f))),
        "dt_first_step": float(dt_s),
    }


def _solve_periodic_steady_exact_newton(L, Nx, params, branch="lp0", shift_half=False, nbar_mode="exact_1"):
    """Exact analytic Newton-PTC solver for RHS == 0 of 1D/PeriodicBoundaries/functions.jl
    on even-symmetric periodic states (used for LP0 plateau states where H-preconditioned
    Newton-Krylov stalls near the saddle-node/Hopf boundary).
    """
    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    dx = L / Nx
    x = -0.5 * L + np.arange(Nx, dtype=float) * dx
    k = 2.0 * math.pi * rfftfreq(Nx, d=dx)
    ik = 1j * k
    k2 = k * k
    vel_op = ik / (1.0 + k2)

    nbar = (Nx / (Nx - 1.0)) if nbar_mode == "trapz_grid" else 1.0
    CH, NaH, NiH = _solve_hss(A_par, Kd, Om0, Om, Omd0, Omd, nbar)

    if shift_half:
        M = Nx // 2
        Mo = M
        E = np.zeros((Nx, M), dtype=float)
        O = np.zeros((Nx, Mo), dtype=float)
        for j in range(M):
            E[Nx // 2 + j, j] = 1.0
            E[Nx // 2 - 1 - j, j] = 1.0
            O[Nx // 2 + j, j] = 1.0
            O[Nx // 2 - 1 - j, j] = -1.0
        e_rows = np.arange(Nx // 2, Nx)
        o_rows = e_rows
        w_e = np.full(M, 2.0 / Nx, dtype=float)
        o_slice = slice(0, M)
    else:
        M = Nx // 2 + 1
        Mo = M - 2
        E = np.zeros((Nx, M), dtype=float)
        E[Nx // 2, 0] = 1.0
        for j in range(1, M - 1):
            E[Nx // 2 + j, j] = 1.0
            E[Nx // 2 - j, j] = 1.0
        E[0, M - 1] = 1.0
        e_rows = np.concatenate([np.arange(Nx // 2, Nx), [0]])
        O = np.zeros((Nx, Mo), dtype=float)
        for j in range(1, M - 1):
            O[Nx // 2 + j, j - 1] = 1.0
            O[Nx // 2 - j, j - 1] = -1.0
        o_rows = np.arange(Nx // 2 + 1, Nx)
        w_e = np.full(M, 2.0 / Nx, dtype=float)
        w_e[0] = 1.0 / Nx
        w_e[-1] = 1.0 / Nx
        o_slice = slice(1, M - 1)

    D2_ee = irfft(rfft(E, axis=0) * (-k2[:, None]), n=Nx, axis=0)[e_rows, :]
    V_oe = irfft(rfft(E, axis=0) * vel_op[:, None], n=Nx, axis=0)[o_rows, :]
    D1_eo = irfft(rfft(O, axis=0) * ik[:, None], n=Nx, axis=0)[e_rows, :]

    xc = -0.5 * dx if shift_half else 0.0
    xe = np.abs(x[e_rows] - xc)
    if branch == "spike":
        C = CH + 0.52 * np.exp(-(xe * xe) / 0.25)
        Na = NaH + 0.12 * np.exp(-(xe * xe) / 0.25)
        Ni = NiH + 0.25 * np.exp(-(xe * xe) / 0.25)
    else:
        bump = 0.5 * (1.0 - np.tanh((xe - 1.5) / 0.3))
        C = CH + 0.70 * bump
        Na = NaH + 0.40 * bump
        Ni = NiH + 1.50 * bump
    Ni += nbar - float(np.dot(w_e, Na + Ni))
    u = np.concatenate([C, Na, Ni])

    inv_dt = 0.5
    eye_M = np.eye(M, dtype=float)
    o_start = o_slice.start
    o_stop = o_slice.stop

    for _ in range(35):
        C, Na, Ni = u[:M], u[M : 2 * M], u[2 * M :]
        sig = Z * C * C - B * C * C * C
        dsig = 2.0 * Z * C - 3.0 * B * C * C
        v = V_oe @ sig
        Co, Nao, Nio = C[o_slice], Na[o_slice], Ni[o_slice]
        ra = -(Omd0 + Omd * C) * Na + Om0 * Ni * (1.0 + Om * Na * Na)
        rC = -D1_eo @ (Co * v) + Dc * (D2_ee @ C) - Kd * C + A_par * Na
        rNa = -D1_eo @ (Nao * v) + Da * (D2_ee @ Na) + ra
        rNi = -D1_eo @ (Nio * v) + (D2_ee @ Ni) - ra
        rNi[-1] = float(np.dot(w_e, Na + Ni)) - nbar
        res = np.concatenate([rC, rNa, rNi])
        nr = float(np.max(np.abs(res)))
        if nr < 1e-11:
            break

        J = np.zeros((3 * M, 3 * M), dtype=float)
        DV_dC = V_oe * dsig[None, :]
        J[:M, :M] = -(D1_eo * Co[None, :]) @ DV_dC + Dc * D2_ee - Kd * eye_M
        J[:M, o_start:o_stop] -= D1_eo * v[None, :]
        J[:M, M : 2 * M] = A_par * eye_M

        dRa_dC = -Omd * Na
        dRa_dNa = -(Omd0 + Omd * C) + 2.0 * Om0 * Om * Ni * Na
        dRa_dNi = Om0 * (1.0 + Om * Na * Na)
        J[M : 2 * M, :M] = -(D1_eo * Nao[None, :]) @ DV_dC + np.diag(dRa_dC)
        J[M : 2 * M, M : 2 * M] = Da * D2_ee + np.diag(dRa_dNa)
        J[M : 2 * M, M + o_start : M + o_stop] -= D1_eo * v[None, :]
        J[M : 2 * M, 2 * M :] = np.diag(dRa_dNi)

        J[2 * M :, :M] = -(D1_eo * Nio[None, :]) @ DV_dC - np.diag(dRa_dC)
        J[2 * M :, M : 2 * M] = -np.diag(dRa_dNa)
        J[2 * M :, 2 * M :] = D2_ee - np.diag(dRa_dNi)
        J[2 * M :, 2 * M + o_start : 2 * M + o_stop] -= D1_eo * v[None, :]

        J[3 * M - 1, :] = 0.0
        J[3 * M - 1, M : 2 * M] = w_e
        J[3 * M - 1, 2 * M :] = w_e
        if inv_dt > 0.0:
            idx = np.arange(3 * M - 1)
            J[idx, idx] -= inv_dt

        step = np.linalg.solve(J, -res)
        alpha = min(1.0, 0.4 / max(0.4, float(np.max(np.abs(step))))) if inv_dt > 0.0 else 1.0
        u = u + alpha * step
        if inv_dt > 0.0:
            if nr < 0.05:
                inv_dt = 0.0
            else:
                inv_dt *= 0.35

    Cf = E @ u[:M]
    Naf = E @ u[M : 2 * M]
    Nif = E @ u[2 * M :]
    F_full = np.vstack([Cf, Naf, Nif])
    rhs_full, _ = _periodic_rhs(F_full, k, ik, k2, vel_op, params)
    rmax = float(np.max(np.abs(rhs_full)))
    return x, Cf, Naf, Nif, CH, rmax


def solve_periodic_steady(L, Nx, params, branch="spike", shift_half=False, nbar_mode="exact_1"):
    """Solve the exact steady-state equation RHS == 0 of 1D/PeriodicBoundaries/functions.jl
    on even-symmetric states.
    Returns (x, C, Na, Ni, CH, rmax).
    """
    if branch != "spike" and Nx <= 1024:
        return _solve_periodic_steady_exact_newton(
            L, Nx, params, branch=branch, shift_half=shift_half, nbar_mode=nbar_mode
        )

    Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B = params
    dx = L / Nx
    x = -0.5 * L + np.arange(Nx, dtype=float) * dx
    k = 2.0 * math.pi * rfftfreq(Nx, d=dx)
    ik = 1j * k
    k2 = k * k
    vel_op = ik / (1.0 + k2)

    nbar = (Nx / (Nx - 1.0)) if nbar_mode == "trapz_grid" else 1.0
    CH, NaH, NiH = _solve_hss(A_par, Kd, Om0, Om, Omd0, Omd, nbar)
    M = (Nx // 2) if shift_half else (Nx // 2 + 1)

    def sym_expand(h):
        out = np.empty(Nx, dtype=float)
        if shift_half:
            out[Nx // 2 :] = h
            out[: Nx // 2] = h[::-1]
        else:
            out[Nx // 2 :] = h[:-1]
            out[0] = h[-1]
            out[1 : Nx // 2] = h[-2:0:-1]
        return out

    def sym_extract(arr):
        if shift_half:
            return arr[Nx // 2 :].copy()
        return np.concatenate([arr[Nx // 2 :], arr[:1]])

    def rhs_3(C, Na, Ni):
        d2C = irfft(rfft(C) * (-k2), n=Nx)
        d2Na = irfft(rfft(Na) * (-k2), n=Nx)
        d2Ni = irfft(rfft(Ni) * (-k2), n=Nx)
        sig = Z * C * C - B * C * C * C
        v = irfft(rfft(sig) * vel_op, n=Nx)
        ac = -irfft(rfft(C * v) * ik, n=Nx)
        aa = -irfft(rfft(Na * v) * ik, n=Nx)
        ai = -irfft(rfft(Ni * v) * ik, n=Nx)
        ra = -(Omd0 + Omd * C) * Na + Om0 * Ni * (1.0 + Om * Na * Na)
        return (
            ac + Dc * d2C - Kd * C + A_par * Na,
            aa + Da * d2Na + ra,
            ai + d2Ni - ra,
        )

    def unpack(u):
        C = sym_expand(u[:M])
        Na = sym_expand(u[M : 2 * M])
        Ni0 = sym_expand(u[2 * M :])
        Ni = Ni0 + (nbar - float(np.mean(Na)) - float(np.mean(Ni0)))
        return C, Na, Ni

    def F_res(u):
        C, Na, Ni = unpack(u)
        rC, rNa, rNi = rhs_3(C, Na, Ni)
        mean_Ni0 = float(np.mean(sym_expand(u[2 * M :])))
        rNi = rNi - float(np.mean(rNi)) - mean_Ni0
        return np.concatenate([sym_extract(rC), sym_extract(rNa), sym_extract(rNi)])

    # Exact linear Fourier preconditioner at H
    dPiH = 2.0 * Z * CH - 3.0 * B * CH * CH
    n_k = len(k)
    fac = (k2 / (1.0 + k2)) * dPiH
    Jblk = np.zeros((n_k, 3, 3), dtype=float)
    Jblk[:, 0, 0] = fac * CH - Dc * k2 - Kd
    Jblk[:, 0, 1] = A_par
    Jblk[:, 1, 0] = fac * NaH - Omd * NaH
    Jblk[:, 1, 1] = -Da * k2 - (Omd0 + Omd * CH) + 2.0 * Om0 * Om * NaH * NiH
    Jblk[:, 1, 2] = Om0 * (1.0 + Om * NaH * NaH)
    Jblk[:, 2, 0] = fac * NiH + Omd * NaH
    Jblk[:, 2, 1] = (Omd0 + Omd * CH) - 2.0 * Om0 * Om * NaH * NiH
    Jblk[:, 2, 2] = -k2 - Om0 * (1.0 + Om * NaH * NaH)
    Jblk[0, 2, :] = [0.0, 0.0, -1.0]
    Jinv = np.linalg.inv(Jblk)

    def precond(v):
        h = np.stack(
            [rfft(sym_expand(v[:M])), rfft(sym_expand(v[M : 2 * M])), rfft(sym_expand(v[2 * M :]))],
            axis=-1,
        )
        out_h = np.einsum("kij,kj->ki", Jinv, h)
        return np.concatenate(
            [
                sym_extract(irfft(out_h[:, 0], n=Nx)),
                sym_extract(irfft(out_h[:, 1], n=Nx)),
                sym_extract(irfft(out_h[:, 2], n=Nx)),
            ]
        )

    M_op = LinearOperator((3 * M, 3 * M), matvec=precond, dtype=float)

    xc = -0.5 * dx if shift_half else 0.0
    C = CH + 0.52 * np.exp(-((x - xc) ** 2) / 0.25)
    Na = NaH + 0.12 * np.exp(-((x - xc) ** 2) / 0.25)
    Ni = NiH + 0.25 * np.exp(-((x - xc) ** 2) / 0.25)
    Ni += nbar - float(np.mean(Na + Ni))
    dt = 0.05
    den_c = 1.0 + dt * (Dc * k2 + Kd)
    den_a = 1.0 + dt * Da * k2
    den_i = 1.0 + dt * k2
    for _ in range(150):
        sig = Z * C * C - B * C * C * C
        v = irfft(rfft(sig) * vel_op, n=Nx)
        ac = -irfft(rfft(C * v) * ik, n=Nx)
        aa = -irfft(rfft(Na * v) * ik, n=Nx)
        ai = -irfft(rfft(Ni * v) * ik, n=Nx)
        ra = -(Omd0 + Omd * C) * Na + Om0 * Ni * (1.0 + Om * Na * Na)
        C = irfft((rfft(C) + dt * rfft(ac + A_par * Na)) / den_c, n=Nx)
        Na = irfft((rfft(Na) + dt * rfft(aa + ra)) / den_a, n=Nx)
        Ni = irfft((rfft(Ni) + dt * rfft(ai - ra)) / den_i, n=Nx)

    u0 = np.concatenate([sym_extract(C), sym_extract(Na), sym_extract(Ni - float(np.mean(Ni)))])
    u = newton_krylov(F_res, u0, f_tol=1e-9, inner_M=M_op, maxiter=40)
    C, Na, Ni = unpack(u)
    rC, rNa, rNi = rhs_3(C, Na, Ni)
    rmax = float(max(np.max(np.abs(rC)), np.max(np.abs(rNa)), np.max(np.abs(rNi))))
    return x, C, Na, Ni, CH, rmax


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
    dt_init = float(inputs["dt_init"])

    solver = str(settings.get("solver", "newton_krylov")).lower()
    branch = str(settings.get("branch", "spike")).lower()
    shift_half = bool(settings.get("shift_half", False))
    nbar_mode = str(settings.get("nbar_mode", "exact_1")).lower()
    n_grid_scale = float(settings.get("n_grid_scale", 1.0))
    dt_imex = float(settings.get("dt_imex", 0.05))
    adapt_steps = int(round(float(settings.get("adapt_steps", 5))))
    seed = int(round(float(settings.get("seed", 123))))

    L = L_base * L_mult
    Nx = int(round(N_base * L_mult * n_grid_scale))
    if Nx % 2 == 1:
        Nx += 1
    params = (Dc, Da, A_par, Kd, Om0, Om, Omd0, Omd, Z, B)

    if solver == "adaptive_midpoint":
        return _run_adaptive_midpoint(L, Nx, dt_init, tol, params, adapt_steps=adapt_steps)
    if solver == "imex_wave":
        ic_mode = "noise" if Z == 0.0 else "bump"
        return _run_imex_time(L, Nx, T_end, dt_imex, params, ic_mode=ic_mode, seed=seed)

    x, C, Na, Ni, CH, rmax = solve_periodic_steady(
        L, Nx, params, branch=branch, shift_half=shift_half, nbar_mode=nbar_mode
    )
    k = 2.0 * math.pi * rfftfreq(Nx, d=L / Nx)
    dC = irfft(rfft(C) * (1j * k), n=Nx)
    return {
        "C_HSS": CH,
        "C_min": float(np.min(C)),
        "C_max": float(np.max(C)),
        "C_center": float(C[Nx // 2]),
        "grad_c_sq_mean": float(np.mean(dC * dC)),
        "mean_N": float(np.mean(Na + Ni)),
        "max_rhs": rmax,
        "dt_first_step": dt_init,
    }
