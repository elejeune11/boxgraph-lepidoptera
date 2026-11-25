import numpy as np


def _etd1_2x2_mats(a11, a12, a21, a22, dt):
    tr = a11 + a22
    half_tr = 0.5 * tr
    diff = 0.5 * (a11 - a22)
    disc = diff * diff + a12 * a21
    s = np.sqrt(disc.astype(np.complex128))
    l1 = half_tr + s
    l2 = half_tr - s
    e1 = np.exp(dt * l1)
    e2 = np.exp(dt * l2)

    def phi1(z):
        big = np.abs(z) > 1e-6
        safe_z = np.where(big, z, 1.0)
        return np.where(big, (np.exp(z) - 1.0) / safe_z, 1.0 + 0.5 * z + z * z / 6.0)

    p1 = dt * phi1(dt * l1)
    p2 = dt * phi1(dt * l2)

    def mat_func(f1, f2, df_mid):
        favg = 0.5 * (f1 + f2)
        bigs = np.abs(s) > 1e-7
        safe_2s = np.where(bigs, 2.0 * s, 1.0)
        fdiff = np.where(bigs, (f1 - f2) / safe_2s, df_mid)
        return (
            (favg + fdiff * diff).real,
            (fdiff * a12).real,
            (fdiff * a21).real,
            (favg - fdiff * diff).real,
        )

    E = mat_func(e1, e2, dt * np.exp(dt * half_tr))
    z0 = dt * half_tr
    bigz = np.abs(z0) > 1e-5
    safe_z0 = np.where(bigz, z0, 1.0)
    dphi = dt * dt * np.where(
        bigz,
        ((z0 - 1.0) * np.exp(z0) + 1.0) / (safe_z0 * safe_z0),
        0.5 + z0 / 3.0 + z0 * z0 / 8.0,
    )
    P = mat_func(p1, p2, dphi)
    return E, P


def run(inputs, settings):
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    chi = float(inputs["chi"])
    chip = float(inputs["chip"])
    alpha = float(inputs["alpha"])
    kappa = float(inputs["kappa"])
    p1b = float(inputs["phi1_bar"])
    p2b = float(inputs["phi2_bar"])

    N = int(settings.get("N", 64))
    L = float(settings.get("L", 2.0))
    dt = float(settings.get("dt", 1e-3))
    T = float(settings.get("T", 6.0))
    t_fit_lo = float(settings.get("t_fit_lo", 0.15))
    t_fit_hi = float(settings.get("t_fit_hi", 2.5))
    S = float(settings.get("S", 0.8))
    init_amp = float(settings.get("init_amp", 0.01))
    seed = int(settings.get("seed", 42))

    rng = np.random.default_rng(seed)
    u1 = init_amp * rng.standard_normal((N, N))
    u2 = init_amp * rng.standard_normal((N, N))
    u1 -= u1.mean()
    u2 -= u2.mean()

    qx = 2.0 * np.pi * np.fft.fftfreq(N, d=L / N)
    qy = 2.0 * np.pi * np.fft.rfftfreq(N, d=L / N)
    QX, QY = np.meshgrid(qx, qy, indexing="ij")
    Q2 = QX * QX + QY * QY
    Q4 = Q2 * Q2
    Q = np.sqrt(Q2)

    kmax = (N // 3) * (2.0 * np.pi / L)
    dealias = (np.abs(QX) <= kmax) & (np.abs(QY) <= kmax)

    # Jacobian of bulk mu^neq at (p1b, p2b), matching Eq. (8)
    j11 = 2.0 * (p1b - c11) ** 2 + 8.0 * (p1b - c11) * (p1b - c12) + 2.0 * (p1b - c12) ** 2 + 2.0 * p2b * p2b * chip
    j22 = 2.0 * (p2b - c21) ** 2 + 8.0 * (p2b - c21) * (p2b - c22) + 2.0 * (p2b - c22) ** 2 + 2.0 * p1b * p1b * chip
    j12 = (chi + alpha) + 4.0 * p1b * p2b * chip
    j21 = (chi - alpha) + 4.0 * p1b * p2b * chip

    # Linear operator L(q) with Eyre stabilization S on the diagonal
    a11 = -(j11 + S) * Q2 - kappa * Q4
    a22 = -(j22 + S) * Q2 - kappa * Q4
    a12 = -j12 * Q2
    a21 = -j21 * Q2
    (e11, e12, e21, e22), (p11, p12, p21, p22) = _etd1_2x2_mats(a11, a12, a21, a22, dt)
    fac = -Q2 * dealias
    p11 *= fac
    p12 *= fac
    p21 *= fac
    p22 *= fac

    # Reference bulk chemical potentials at (p1b, p2b)
    m1_0 = 2.0 * (p1b - c11) * (p1b - c12) * (2.0 * p1b - c11 - c12) + (chi + alpha) * p2b + 2.0 * chip * p1b * p2b * p2b
    m2_0 = 2.0 * (p2b - c21) * (p2b - c22) * (2.0 * p2b - c21 - c22) + (chi - alpha) * p1b + 2.0 * chip * p2b * p1b * p1b

    u1h = np.fft.rfft2(u1)
    u2h = np.fft.rfft2(u2)
    nsteps = max(1, int(round(T / dt)))

    w = np.full_like(Q, 2.0)
    w[:, 0] = 1.0
    if N % 2 == 0:
        w[:, -1] = 1.0
    mask = Q > 0.0

    t_start = min(0.05, 0.5 * t_fit_lo)
    check_steps = set(int(round(t / dt)) for t in np.geomspace(t_start, T, 32)) | {nsteps}
    ts = []
    Rs = []
    J_mag_last = 0.0
    J_comp_last = 0.0

    for step in range(1, nsteps + 1):
        p1 = p1b + u1
        p2 = p2b + u2
        m1 = 2.0 * (p1 - c11) * (p1 - c12) * (2.0 * p1 - c11 - c12) + (chi + alpha) * p2 + 2.0 * chip * p1 * p2 * p2
        m2 = 2.0 * (p2 - c21) * (p2 - c22) * (2.0 * p2 - c21 - c22) + (chi - alpha) * p1 + 2.0 * chip * p2 * p1 * p1
        # Nonlinear residual after subtracting linearized part + Eyre S
        g1 = (m1 - m1_0) - (j11 + S) * u1 - j12 * u2
        g2 = (m2 - m2_0) - j21 * u1 - (j22 + S) * u2
        g1h = np.fft.rfft2(g1)
        g2h = np.fft.rfft2(g2)
        u1h_new = e11 * u1h + e12 * u2h + p11 * g1h + p12 * g2h
        u2h_new = e21 * u1h + e22 * u2h + p21 * g1h + p22 * g2h
        u1h, u2h = u1h_new, u2h_new
        u1 = np.fft.irfft2(u1h, s=(N, N))
        u2 = np.fft.irfft2(u2h, s=(N, N))

        if step in check_steps:
            S2 = np.abs(u1h) ** 2
            denom_R = float(np.sum(w[mask] * S2[mask]))
            R_val = float(2.0 * np.pi * np.sum(w[mask] * S2[mask] / Q[mask]) / (denom_R + 1e-30))
            ts.append(step * dt)
            Rs.append(R_val)
            if step == nsteps:
                du1x = np.fft.irfft2(1j * QX * u1h, s=(N, N))
                du1y = np.fft.irfft2(1j * QY * u1h, s=(N, N))
                du2x = np.fft.irfft2(1j * QX * u2h, s=(N, N))
                du2y = np.fft.irfft2(1j * QY * u2h, s=(N, N))
                Jx = float(np.mean(u1 * du2x - u2 * du1x))
                Jy = float(np.mean(u1 * du2y - u2 * du1y))
                J_mag_last = float(np.hypot(Jx, Jy))
                J_comp_last = 0.5 * (abs(Jx) + abs(Jy))

    ts_arr = np.array(ts)
    Rs_arr = np.array(Rs)
    fit_mask = (ts_arr >= t_fit_lo) & (ts_arr <= t_fit_hi) & (Rs_arr > 0.0)
    if np.sum(fit_mask) >= 3:
        beta = float(np.polyfit(np.log(ts_arr[fit_mask]), np.log(Rs_arr[fit_mask]), 1)[0])
    else:
        beta = float("nan")

    return {
        "coarsening_exponent": beta,
        "R_final": float(Rs_arr[-1]),
        "J_mag": J_mag_last,
        "J_comp_avg": J_comp_last,
    }
