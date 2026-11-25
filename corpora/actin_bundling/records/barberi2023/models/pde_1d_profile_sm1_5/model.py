import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq


def _solve_H(A, Kd, Om0, Om, Omd0, Omd, Nbar):
    ratio = A / Kd
    na = float(
        brentq(
            lambda x: Om0 * (1.0 + Om * x * x) * (Nbar - x) - (Omd0 + Omd * ratio * x) * x,
            0.0,
            Nbar,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    return ratio * na, na, Nbar - na


def _tri(diag_val, off_val, n):
    ab = np.zeros((3, n), dtype=float)
    ab[0, 1:] = off_val
    ab[1, :] = diag_val
    ab[2, :-1] = off_val
    return ab


def _eval_res_fn(L, N_full, Dc, Da, Di, A, Kd, Om0, Om, Omd0, Omd, Z, B, Nbar):
    M = N_full // 2
    dx = L / N_full
    inv_dx = 1.0 / dx
    inv_dx2 = inv_dx * inv_dx
    nv = 5

    def eval_res(y):
        Y = y.reshape(M, nv)
        C, Na, Ni, V, Q = Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3], Y[:, 4]

        def divJ(S, D):
            J = np.zeros(M + 1, dtype=float)
            J[1:M] = 0.5 * V[:-1] * (S[:-1] + S[1:]) - D * (S[1:] - S[:-1]) * inv_dx
            return (J[1:] - J[:-1]) * inv_dx

        R = Om0 * (1.0 + Om * Na * Na) * Ni - (Omd0 + Omd * C) * Na
        rC = -divJ(C, Dc) + A * Na - Kd * C
        rNa = -divJ(Na, Da) + R
        rNi = -divJ(Ni, Di) - R
        # Replace redundant boundary flux balance with exact global nucleator conservation
        rNi[-1] = Q[-1] - Nbar

        Pi = Z * C * C - B * C * C * C
        rV = np.zeros(M, dtype=float)
        V_left = np.empty(M - 1, dtype=float)
        V_left[0] = 0.0
        V_left[1:] = V[:-2]
        rV[:-1] = V[:-1] - (V[1:] - 2.0 * V[:-1] + V_left) * inv_dx2 - (Pi[1:] - Pi[:-1]) * inv_dx
        rV[-1] = V[-1]

        rQ = np.empty(M, dtype=float)
        rQ[0] = Q[0] - (Na[0] + Ni[0]) / M
        rQ[1:] = Q[1:] - Q[:-1] - (Na[1:] + Ni[1:]) / M

        out = np.empty((M, nv), dtype=float)
        out[:, 0] = rC
        out[:, 1] = rNa
        out[:, 2] = rNi
        out[:, 3] = rV
        out[:, 4] = rQ
        return out.ravel()

    return eval_res


def _build_banded_jac(eval_res, y, r0, M):
    nv = 5
    bw = 9
    Ntot = nv * M
    ab = np.zeros((2 * bw + 1, Ntot), dtype=float)
    eps = 1e-7
    for v in range(nv):
        for mod in range(3):
            cols = np.arange(mod, M, 3) * nv + v
            dy = np.zeros(Ntot, dtype=float)
            dy[cols] = eps
            dr = (eval_res(y + dy) - r0) / eps
            for shift in range(-bw, bw + 1):
                rows = cols + shift
                valid = (rows >= 0) & (rows < Ntot)
                if np.any(valid):
                    ab[bw + shift, cols[valid]] = dr[rows[valid]]
    return ab


def _pack_state(L, N_full, C, Na, Ni, Z, B):
    M = N_full // 2
    dx = L / N_full
    inv_dx = 1.0 / dx
    inv_dx2 = inv_dx * inv_dx
    Pi = Z * C * C - B * C * C * C
    abV = _tri(1.0 + 2.0 * inv_dx2, -inv_dx2, M - 1)
    V = np.zeros(M, dtype=float)
    V[:-1] = solve_banded((1, 1), abV, (Pi[1:] - Pi[:-1]) * inv_dx)
    Q = np.cumsum(Na + Ni) / M
    return np.column_stack([C, Na, Ni, V, Q]).ravel()


def _newton_ptc(L, N_full, y, params, inv_dt0=0.0, maxit=30):
    Dc, Da, Di, A, Kd, Om0, Om, Omd0, Omd, Z, B, Nbar = params
    M = N_full // 2
    nv = 5
    bw = 9
    eval_res = _eval_res_fn(L, N_full, Dc, Da, Di, A, Kd, Om0, Om, Omd0, Omd, Z, B, Nbar)
    c_idx = np.arange(M) * nv + 0
    na_idx = np.arange(M) * nv + 1
    ni_idx = np.arange(M - 1) * nv + 2

    inv_dt = inv_dt0
    nr = float("inf")
    for _ in range(maxit):
        r0 = eval_res(y)
        nr = float(np.max(np.abs(r0)))
        if nr < 1e-11:
            return y, nr
        ab = _build_banded_jac(eval_res, y, r0, M)
        if inv_dt > 0.0:
            ab[bw, c_idx] -= inv_dt
            ab[bw, na_idx] -= inv_dt
            ab[bw, ni_idx] -= inv_dt
        step = solve_banded((bw, bw), ab, -r0)
        if inv_dt > 0.0:
            max_s = float(np.max(np.abs(step)))
            alpha = min(1.0, 0.3 / max(0.3, max_s))
            yn = y + alpha * step
            nr_n = float(np.max(np.abs(eval_res(yn))))
            y = yn
            if nr_n < 0.05:
                inv_dt = 0.0
            elif nr_n < nr:
                inv_dt *= 0.35
            else:
                inv_dt = max(inv_dt * 0.8, 1e-3)
        else:
            n0 = float(np.linalg.norm(r0))
            alpha = 1.0
            for _ in range(10):
                yn = y + alpha * step
                if float(np.linalg.norm(eval_res(yn))) < n0 or alpha < 1e-2:
                    y = yn
                    break
                alpha *= 0.5
    return y, float(np.max(np.abs(eval_res(y))))


def _relax_imex(L, N_sim, C0, Na0, Ni0, params, T_end=220.0, dt=0.025):
    Dc, Da, Di, A, Kd, Om0, Om, Omd0, Omd, Z, B, Nbar = params
    M = N_sim // 2
    dx = L / N_sim
    idx = 1.0 / dx
    idx2 = idx * idx
    abV = _tri(1.0 + 2.0 * idx2, -idx2, M - 1)

    def mdiff(D, kl):
        ab = np.zeros((3, M), dtype=float)
        ab[0, 1:] = -dt * D * idx2
        ab[2, :-1] = -dt * D * idx2
        d = np.full(M, 1.0 + dt * kl + 2.0 * dt * D * idx2, dtype=float)
        d[0] -= dt * D * idx2
        d[-1] -= dt * D * idx2
        ab[1, :] = d
        return ab

    abC = mdiff(Dc, Kd)
    abNa = mdiff(Da, 0.0)
    abNi = mdiff(Di, 0.0)
    C, Na, Ni = C0.copy(), Na0.copy(), Ni0.copy()
    nsteps = int(round(T_end / dt))
    for s in range(nsteps):
        Pi = Z * C * C - B * C * C * C
        Vf = solve_banded(
            (1, 1), abV, (Pi[1:] - Pi[:-1]) * idx, overwrite_ab=False, overwrite_b=True, check_finite=False
        )
        Jc = 0.5 * Vf * (C[:-1] + C[1:])
        Jna = 0.5 * Vf * (Na[:-1] + Na[1:])
        Jni = 0.5 * Vf * (Ni[:-1] + Ni[1:])
        R = Om0 * (1.0 + Om * Na * Na) * Ni - (Omd0 + Omd * C) * Na
        rC = C + dt * A * Na
        rC[:-1] -= dt * Jc * idx
        rC[1:] += dt * Jc * idx
        rNa = Na + dt * R
        rNa[:-1] -= dt * Jna * idx
        rNa[1:] += dt * Jna * idx
        rNi = Ni - dt * R
        rNi[:-1] -= dt * Jni * idx
        rNi[1:] += dt * Jni * idx
        Cn = solve_banded((1, 1), abC, rC, overwrite_ab=False, overwrite_b=True, check_finite=False)
        Nan = solve_banded((1, 1), abNa, rNa, overwrite_ab=False, overwrite_b=True, check_finite=False)
        Nin = solve_banded((1, 1), abNi, rNi, overwrite_ab=False, overwrite_b=True, check_finite=False)
        if s % 300 == 299 and float(np.max(np.abs(Cn - C))) / dt < 3e-4:
            C, Na, Ni = Cn, Nan, Nin
            break
        C, Na, Ni = Cn, Nan, Nin
    return C, Na, Ni


def _continue_omega(L, N_full, y, om_start, om_target, base_params):
    Dc, Da, Di, A, Kd, Om0, _, Omd0, Omd, Z, B, Nbar = base_params
    if abs(om_target - om_start) < 1e-12:
        return _newton_ptc(L, N_full, y, base_params, inv_dt0=0.0, maxit=20)
    n_steps = max(2, int(math.ceil(abs(om_target - om_start) / 0.35)))
    for om_s in np.linspace(om_start, om_target, n_steps + 1)[1:]:
        p = (Dc, Da, Di, A, Kd, Om0, float(om_s), Omd0, Omd, Z, B, Nbar)
        y, nr = _newton_ptc(L, N_full, y, p, inv_dt0=0.0, maxit=20)
    return y, nr



def _solve_branch_full(L, N_full, branch, params):
    Dc, Da, Di, A, Kd, Om0, Om, Omd0, Omd, Z, B, Nbar = params
    if N_full % 2 == 1:
        N_full += 1
    M = N_full // 2
    dx = L / N_full
    x = (np.arange(M, dtype=float) + 0.5) * dx
    c_h, na_h, ni_h = _solve_H(A, Kd, Om0, Om, Omd0, Omd, Nbar)

    if branch == "h":
        C = np.full(M, c_h)
        Na = np.full(M, na_h)
        Ni = np.full(M, ni_h)
        V = np.zeros(M)
        return x, dx, C, Na, Ni, V, 0.0

    if branch == "spike":
        C0 = c_h + 0.52 * np.exp(-x * x / 0.5)
        Na0 = na_h + 0.12 * np.exp(-x * x / 0.5)
        Ni0 = ni_h + 0.25 * np.exp(-x * x / 0.5)
        Ni0 += Nbar - float(np.mean(Na0 + Ni0))
        y0 = _pack_state(L, N_full, C0, Na0, Ni0, Z, B)
        y, nr = _newton_ptc(L, N_full, y0, params, inv_dt0=0.5, maxit=35)
    elif branch == "lp0":
        bump = 0.5 * (1.0 - np.tanh((x - 1.5) / 0.3))
        C0 = c_h + 0.70 * bump
        Na0 = na_h + 0.40 * bump
        Ni0 = ni_h + 1.50 * bump
        Ni0 += Nbar - float(np.mean(Na0 + Ni0))
        y0 = _pack_state(L, N_full, C0, Na0, Ni0, Z, B)
        y, nr = _newton_ptc(L, N_full, y0, params, inv_dt0=1.0, maxit=40)
    else:
        N_sim = 256
        Ms = N_sim // 2
        xs = (np.arange(Ms, dtype=float) + 0.5) * (L / N_sim)

        def plateau_seed(hw):
            env = 0.5 * (1.0 - np.tanh((xs - hw) / 0.3))
            C_s = 0.22 + 0.65 * env
            Na_s = 0.22 + 0.60 * env
            Ni_s = 0.78 + 0.80 * env
            scale = Nbar / float(np.mean(Na_s + Ni_s))
            return C_s, Na_s * scale, Ni_s * scale

        def to_fine_and_polish(Cs, Nas, Nis, om_s):
            Cf = np.interp(x, xs, Cs)
            Naf = np.interp(x, xs, Nas)
            Nif = np.interp(x, xs, Nis)
            Nif += Nbar - float(np.mean(Naf + Nif))
            yf0 = _pack_state(L, N_full, Cf, Naf, Nif, Z, B)
            p_s = (Dc, Da, Di, A, Kd, Om0, om_s, Omd0, Omd, Z, B, Nbar)
            yf, _ = _newton_ptc(L, N_full, yf0, p_s, inv_dt0=0.0, maxit=25)
            return _continue_omega(L, N_full, yf, om_s, Om, params)

        if branch == "lp2":
            om_s = 10.0
            p_s = (Dc, Da, Di, A, Kd, Om0, om_s, Omd0, Omd, Z, B, Nbar)
            Cs, Nas, Nis = _relax_imex(L, N_sim, *plateau_seed(3.5), p_s)
            y, nr = to_fine_and_polish(Cs, Nas, Nis, om_s)
        elif branch == "lp4":
            om_s = 14.0
            p_s = (Dc, Da, Di, A, Kd, Om0, om_s, Omd0, Omd, Z, B, Nbar)
            Cs, Nas, Nis = _relax_imex(L, N_sim, *plateau_seed(5.5), p_s)
            y, nr = to_fine_and_polish(Cs, Nas, Nis, om_s)
        elif branch == "lp6":
            om_s = 15.5
            p_s = (Dc, Da, Di, A, Kd, Om0, om_s, Omd0, Omd, Z, B, Nbar)
            Cs, Nas, Nis = _relax_imex(L, N_sim, *plateau_seed(6.2), p_s)
            y, nr = to_fine_and_polish(Cs, Nas, Nis, om_s)
        elif branch == "lp3":
            p_ref = (Dc, Da, Di, A, Kd, Om0, 14.0, Omd0, Omd, Z, B, Nbar)
            C4, Na4, Ni4 = _relax_imex(L, N_sim, *plateau_seed(5.5), p_ref)
            C3 = np.interp(xs + 1.15, xs, C4, right=C4[-1])
            Na3 = np.interp(xs + 1.15, xs, Na4, right=Na4[-1])
            Ni3 = np.interp(xs + 1.15, xs, Ni4, right=Ni4[-1])
            Ni3 += Nbar - float(np.mean(Na3 + Ni3))
            om_s = 11.0
            p_s = (Dc, Da, Di, A, Kd, Om0, om_s, Omd0, Omd, Z, B, Nbar)
            C3, Na3, Ni3 = _relax_imex(L, N_sim, C3, Na3, Ni3, p_s, T_end=180.0)
            y, nr = to_fine_and_polish(C3, Na3, Ni3, om_s)
        elif branch == "lp5":
            p_ref = (Dc, Da, Di, A, Kd, Om0, 15.5, Omd0, Omd, Z, B, Nbar)
            C6, Na6, Ni6 = _relax_imex(L, N_sim, *plateau_seed(6.2), p_ref)
            C5 = np.interp(xs + 0.874, xs, C6, right=C6[-1])
            Na5 = np.interp(xs + 0.874, xs, Na6, right=Na6[-1])
            Ni5 = np.interp(xs + 0.874, xs, Ni6, right=Ni6[-1])
            Ni5 += Nbar - float(np.mean(Na5 + Ni5))
            om_s = 15.0
            p_s = (Dc, Da, Di, A, Kd, Om0, om_s, Omd0, Omd, Z, B, Nbar)
            C5, Na5, Ni5 = _relax_imex(L, N_sim, C5, Na5, Ni5, p_s, T_end=180.0)
            y, nr = to_fine_and_polish(C5, Na5, Ni5, om_s)
        elif branch == "lp7":
            p_ref = (Dc, Da, Di, A, Kd, Om0, 15.5, Omd0, Omd, Z, B, Nbar)
            C6, Na6, Ni6 = _relax_imex(L, N_sim, *plateau_seed(6.2), p_ref)
            x_src = np.abs(xs - 0.874)
            C7, Na7, Ni7 = np.interp(x_src, xs, C6), np.interp(x_src, xs, Na6), np.interp(x_src, xs, Ni6)
            Ni7 += Nbar - float(np.mean(Na7 + Ni7))
            om_s = 16.5
            p_s = (Dc, Da, Di, A, Kd, Om0, om_s, Omd0, Omd, Z, B, Nbar)
            C7, Na7, Ni7 = _relax_imex(L, N_sim, C7, Na7, Ni7, p_s, T_end=180.0)
            y, nr = to_fine_and_polish(C7, Na7, Ni7, om_s)
        else:
            raise ValueError(f"Unknown branch: {branch!r}")

    Y = y.reshape(M, 5)
    return x, dx, Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3], float(nr)


def _flux_terms(S, V, D, dx):
    M = len(S)
    inv_dx = 1.0 / dx
    Ja = np.zeros(M + 1, dtype=float)
    Jd = np.zeros(M + 1, dtype=float)
    Ja[1:M] = 0.5 * V[:-1] * (S[:-1] + S[1:])
    Jd[1:M] = -D * (S[1:] - S[:-1]) * inv_dx
    neg_dJa = -(Ja[1:] - Ja[:-1]) * inv_dx
    neg_dJd = -(Jd[1:] - Jd[:-1]) * inv_dx
    return neg_dJa, neg_dJd

def run(inputs, settings):
    L = float(inputs["L"])
    params = (
        float(inputs["Dc"]), float(inputs["Da"]), float(inputs["Di"]),
        float(inputs["A"]), float(inputs["Kd"]), float(inputs["Omega0"]),
        float(inputs["Omega"]), float(inputs["Omega_d0"]), float(inputs["Omega_d"]),
        float(inputs["Z"]), float(inputs["B"]), float(inputs["Nbar"]),
    )
    branch = str(settings.get("branch", "lp0")).lower()
    n_grid = int(round(settings.get("n_grid", 1024)))
    align_peak = bool(settings.get("align_peak", False))
    x_max = settings.get("x_max", None)
    abs_x_min = settings.get("abs_x_min", None)

    ref = inputs["ref_C"]
    xs_all = [float(v) for v in ref["x"]]
    ys_all = [float(v) for v in ref["y"]]

    x_center = 0.0
    if align_peak and len(xs_all) >= 3:
        # Align symmetric X=0 peak with the parabolic peak of ref_C (periodic translation invariance)
        imax = int(np.argmax(ys_all))
        if 0 < imax < len(xs_all) - 1:
            x0, x1, x2 = xs_all[imax - 1], xs_all[imax], xs_all[imax + 1]
            y0, y1, y2 = ys_all[imax - 1], ys_all[imax], ys_all[imax + 1]
            c = np.polyfit([x0, x1, x2], [y0, y1, y2], 2)
            x_center = float(-c[1] / (2.0 * c[0])) if abs(c[0]) > 1e-12 else x1
        else:
            x_center = xs_all[imax]

    xs_out = []
    for xv in xs_all:
        if x_max is not None and xv > float(x_max):
            continue
        if abs_x_min is not None and abs(xv) < float(abs_x_min):
            continue
        xs_out.append(xv)

    x, dx, C, Na, Ni, V, _ = _solve_branch_full(L, n_grid, branch, params)
    # Cell-centered V (odd function of X)
    V_left = np.empty_like(V)
    V_left[0] = 0.0
    V_left[1:] = V[:-1]
    V_cell = 0.5 * (V_left + V)

    Jca, Jcd = _flux_terms(C, V, params[0], dx)
    Jaa, Jad = _flux_terms(Na, V, params[1], dx)
    Jia, Jid = _flux_terms(Ni, V, params[2], dx)

    xf = np.concatenate([-x[::-1], x])
    def _even(arr):
        return np.concatenate([arr[::-1], arr])
    def _odd(arr):
        return np.concatenate([-arr[::-1], arr])

    xq = np.array(xs_out, dtype=float) - x_center
    # Periodic wrap into [-L/2, L/2]
    xq = ((xq + 0.5 * L) % L) - 0.5 * L

    def _eval(arr_f):
        return {"x": list(xs_out), "y": [float(v) for v in np.interp(xq, xf, arr_f)]}

    return {
        "C_curve": _eval(_even(C)),
        "Na_curve": _eval(_even(Na)),
        "Ni_curve": _eval(_even(Ni)),
        "V_curve": _eval(_odd(V_cell)),
        "neg_dX_Jca_curve": _eval(_even(Jca)),
        "neg_dX_Jcd_curve": _eval(_even(Jcd)),
        "neg_dX_Jaa_curve": _eval(_even(Jaa)),
        "neg_dX_Jad_curve": _eval(_even(Jad)),
        "neg_dX_Jia_curve": _eval(_even(Jia)),
        "neg_dX_Jid_curve": _eval(_even(Jid)),
    }
