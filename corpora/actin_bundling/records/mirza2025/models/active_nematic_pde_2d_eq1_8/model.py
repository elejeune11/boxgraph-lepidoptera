import math
import numpy as np
from scipy.optimize import minimize_scalar


def run(inputs, settings):
    kd = float(inputs["k_d"])
    Lf = float(inputs["L"])
    c0 = float(inputs["c_0"])
    b = float(inputs["b"])
    er = float(inputs["eta_rot"])
    beta = float(inputs["beta"])
    lo = float(inputs["lambda_o"])
    kap = float(inputs["kappa"])
    ratio = float(inputs["activity_ratio"])

    two_a = c0 + lo
    G = c0 if c0 >= 0.0 else -2.0 * c0
    q0 = 0.0 if c0 >= 0.0 else -0.5 * math.sqrt(-c0 / b)

    def lam_fn(x):
        if x <= 0.0:
            return 1e30
        num = (1.0 / x + 4.0) * (x + kd) * (2.0 * G + 2.0 * Lf * x)
        den = (1.0 + kap * q0) * (2.0 * G + 2.0 * Lf * x) + kap * (beta * (x + kd) + 2.0 * lo * q0)
        return num / den if den > 0.0 else 1e30

    xs = np.geomspace(1e-4, 1e3, 2000)
    vals = np.array([lam_fn(float(x)) for x in xs])
    i0 = int(np.argmin(vals))
    res = minimize_scalar(
        lam_fn,
        bounds=(float(xs[max(0, i0 - 3)]), float(xs[min(len(xs) - 1, i0 + 3)])),
        method="bounded",
        options={"xatol": 1e-12},
    )
    lam_c = float(res.fun)
    nu_c = math.sqrt(float(res.x))
    l_pat = 2.0 * math.pi / nu_c
    lam = ratio * lam_c

    N = int(round(float(settings.get("N_grid", 32))))
    dt = float(settings.get("dt", 0.04))
    t_max = float(settings.get("t_max", 75.0))
    L_mult = float(settings.get("L_mult", 2.0))
    maxiter = int(settings.get("pcg_maxiter", 4))
    seed = int(settings.get("seed", 1))

    Ldom = L_mult * l_pat
    dx = Ldom / N
    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    ky = 2.0 * np.pi * np.fft.rfftfreq(N, d=dx)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K2 = KX * KX + KY * KY
    kmax = np.pi / dx
    dealias = (np.abs(KX) < (2.0 / 3.0) * kmax) & (np.abs(KY) < (2.0 / 3.0) * kmax)

    rfft = np.fft.rfft2
    irfft = lambda a: np.fft.irfft2(a, s=(N, N))

    mu_s = 2.0 - (beta * beta) / er
    a11 = 1.0 + 0.5 * mu_s * K2 + 3.0 * KX * KX
    a12 = 3.0 * KX * KY
    a22 = 1.0 + 0.5 * mu_s * K2 + 3.0 * KY * KY
    det = a11 * a22 - a12 * a12
    inv11 = a22 / det
    inv12 = -a12 / det
    inv22 = a11 / det

    def apply_A(rho_arr, ux, uy):
        uh, vh = rfft(ux), rfft(uy)
        uxx = irfft(1j * KX * uh)
        uxy = irfft(1j * KY * uh)
        uyx = irfft(1j * KX * vh)
        uyy = irfft(1j * KY * vh)
        trd = uxx + uyy
        d11 = 0.5 * (uxx - uyy)
        d12 = 0.5 * (uxy + uyx)
        txx = rho_arr * (mu_s * d11 + 3.0 * trd)
        tyy = rho_arr * (-mu_s * d11 + 3.0 * trd)
        txy = rho_arr * (mu_s * d12)
        txxh, txyh, tyyh = rfft(txx), rfft(txy), rfft(tyy)
        return (
            rho_arr * ux - irfft(1j * KX * txxh + 1j * KY * txyh),
            rho_arr * uy - irfft(1j * KX * txyh + 1j * KY * tyyh),
        )

    def precond(rx, ry, rm):
        rxh, ryh = rfft(rx / rm), rfft(ry / rm)
        return irfft(inv11 * rxh + inv12 * ryh), irfft(inv12 * rxh + inv22 * ryh)

    def solve_v_pcg(rho_arr, bx, by, vx_in, vy_in):
        rm = float(np.mean(rho_arr))
        Ax, Ay = apply_A(rho_arr, vx_in, vy_in)
        rx, ry = bx - Ax, by - Ay
        bnorm = math.sqrt(float(np.sum(bx * bx + by * by))) + 1e-30
        if math.sqrt(float(np.sum(rx * rx + ry * ry))) < 1e-3 * bnorm:
            return vx_in, vy_in
        zx, zy = precond(rx, ry, rm)
        px, py = zx.copy(), zy.copy()
        rz = float(np.sum(rx * zx + ry * zy))
        vx_c, vy_c = vx_in, vy_in
        for _ in range(maxiter):
            Apx, Apy = apply_A(rho_arr, px, py)
            pAp = float(np.sum(px * Apx + py * Apy))
            if pAp <= 0.0:
                break
            al = rz / pAp
            vx_c = vx_c + al * px
            vy_c = vy_c + al * py
            rx = rx - al * Apx
            ry = ry - al * Apy
            if math.sqrt(float(np.sum(rx * rx + ry * ry))) < 1e-3 * bnorm:
                break
            zx, zy = precond(rx, ry, rm)
            rz2 = float(np.sum(rx * zx + ry * zy))
            be = rz2 / rz
            px = zx + be * px
            py = zy + be * py
            rz = rz2
        return vx_c, vy_c

    # Smooth low-mode deterministic initial condition invariant to spatial grid resolution N
    x = np.linspace(0.0, Ldom, N, endpoint=False)
    y = np.linspace(0.0, Ldom, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    rng = np.random.default_rng(seed)
    rho = np.ones((N, N))
    q1 = np.full((N, N), q0)
    q2 = np.zeros((N, N))
    for mx in range(-3, 4):
        for my in range(-3, 4):
            if mx == 0 and my == 0:
                continue
            ph1, ph2, ph3 = rng.uniform(0.0, 2.0 * np.pi, 3)
            amp = 0.015 / (1.0 + mx * mx + my * my)
            arg = 2.0 * np.pi * (mx * X + my * Y) / Ldom
            rho += amp * np.cos(arg + ph1)
            q1 += 0.5 * amp * np.cos(arg + ph2)
            q2 += 0.5 * amp * np.cos(arg + ph3)

    vx = np.zeros((N, N))
    vy = np.zeros((N, N))
    denom_rho = 1.0 + dt * (K2 + kd)
    denom_q = 1.0 + dt * (Lf / er) * K2

    for _ in range(int(round(t_max / dt))):
        rh = rfft(rho)
        q1h = rfft(q1)
        q2h = rfft(q2)
        rx, ry = irfft(1j * KX * rh), irfft(1j * KY * rh)
        q1x, q1y, q1lap = irfft(1j * KX * q1h), irfft(1j * KY * q1h), irfft(-K2 * q1h)
        q2x, q2y, q2lap = irfft(1j * KX * q2h), irfft(1j * KY * q2h), irfft(-K2 * q2h)
        S2 = 4.0 * (q1 * q1 + q2 * q2)
        H1 = (two_a - rho * lo + b * S2) * q1 - Lf * q1lap - (Lf / rho) * (rx * q1x + ry * q1y)
        H2 = (two_a - rho * lo + b * S2) * q2 - Lf * q2lap - (Lf / rho) * (rx * q2x + ry * q2y)
        div_rq1 = rho * q1lap + rx * q1x + ry * q1y
        div_rq2 = rho * q2lap + rx * q2x + ry * q2y
        anti = 2.0 * Lf * (q1 * div_rq2 - q2 * div_rq1)
        sxx_r = rho * (lam * (1.0 + kap * q1) - (beta / er) * H1 - 2.0 * Lf * (q1x * q1x + q2x * q2x))
        syy_r = rho * (lam * (1.0 - kap * q1) + (beta / er) * H1 - 2.0 * Lf * (q1y * q1y + q2y * q2y))
        sxy_r = rho * (lam * kap * q2 - (beta / er) * H2 - 2.0 * Lf * (q1x * q1y + q2x * q2y)) + anti
        syx_r = rho * (lam * kap * q2 - (beta / er) * H2 - 2.0 * Lf * (q1x * q1y + q2x * q2y)) - anti
        fx_r = irfft(1j * KX * rfft(sxx_r) + 1j * KY * rfft(sxy_r))
        fy_r = irfft(1j * KX * rfft(syx_r) + 1j * KY * rfft(syy_r))
        vx, vy = solve_v_pcg(rho, fx_r, fy_r, vx, vy)

        vhx, vhy = rfft(vx), rfft(vy)
        vxx, vxy = irfft(1j * KX * vhx), irfft(1j * KY * vhx)
        vyx, vyy = irfft(1j * KX * vhy), irfft(1j * KY * vhy)
        d11 = 0.5 * (vxx - vyy)
        d12 = 0.5 * (vxy + vyx)
        w12 = 0.5 * (vxy - vyx)
        nl_rho = -irfft(1j * KX * rfft(rho * vx) + 1j * KY * rfft(rho * vy)) + kd
        rho = np.maximum(irfft(rfft(rho + dt * nl_rho) * dealias / denom_rho), 1e-2)

        rhs1 = -(vx * q1x + vy * q1y) + 2.0 * w12 * q2 - (beta / er) * d11 + (1.0 / er) * (
            rho * lo * q1 + (Lf / rho) * (rx * q1x + ry * q1y)
        )
        rhs2 = -(vx * q2x + vy * q2y) - 2.0 * w12 * q1 - (beta / er) * d12 + (1.0 / er) * (
            rho * lo * q2 + (Lf / rho) * (rx * q2x + ry * q2y)
        )
        loc_denom = 1.0 + dt * (two_a + b * S2) / er
        q1 = irfft(rfft((q1 + dt * rhs1) / loc_denom) * dealias / denom_q)
        q2 = irfft(rfft((q2 + dt * rhs2) / loc_denom) * dealias / denom_q)

    rh = rfft(rho)
    rx, ry = irfft(1j * KX * rh), irfft(1j * KY * rh)
    omega = float(np.mean(q1 * (rx * rx - ry * ry) + 2.0 * q2 * rx * ry))
    S = 2.0 * np.sqrt(q1 * q1 + q2 * q2)
    return {
        "omega": omega,
        "rho_max": float(np.max(rho)),
        "rho_min": float(np.min(rho)),
        "S_max": float(np.max(S)),
    }
