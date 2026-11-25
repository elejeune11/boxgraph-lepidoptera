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
    ts_target = [float(v) for v in inputs["target_J"]["x"]]

    N = int(settings.get("N", 48))
    L = float(settings.get("L", 2.0))
    dt = float(settings.get("dt", 4e-3))
    S = float(settings.get("S", 1.5))
    init_amp = float(settings.get("init_amp", 0.01))
    n_seeds = int(settings.get("n_seeds", 3))
    seed0 = int(settings.get("seed", 42))

    T = max(ts_target)
    nsteps = max(1, int(round(T / dt)))
    step_map = {}
    for idx, tv in enumerate(ts_target):
        st = max(1, min(nsteps, int(round(tv / dt))))
        step_map.setdefault(st, []).append(idx)

    qx = 2.0 * np.pi * np.fft.fftfreq(N, d=L / N)
    qy = 2.0 * np.pi * np.fft.rfftfreq(N, d=L / N)
    QX, QY = np.meshgrid(qx, qy, indexing="ij")
    Q2 = QX * QX + QY * QY
    Q4 = Q2 * Q2
    kmax = (N // 3) * (2.0 * np.pi / L)
    dealias = (np.abs(QX) <= kmax) & (np.abs(QY) <= kmax)

    j11 = 2.0 * (p1b - c11) ** 2 + 8.0 * (p1b - c11) * (p1b - c12) + 2.0 * (p1b - c12) ** 2 + 2.0 * p2b * p2b * chip
    j22 = 2.0 * (p2b - c21) ** 2 + 8.0 * (p2b - c21) * (p2b - c22) + 2.0 * (p2b - c22) ** 2 + 2.0 * p1b * p1b * chip
    j12 = (chi + alpha) + 4.0 * p1b * p2b * chip
    j21 = (chi - alpha) + 4.0 * p1b * p2b * chip

    (e11, e12, e21, e22), (p11, p12, p21, p22) = _etd1_2x2_mats(
        -(j11 + S) * Q2 - kappa * Q4, -j12 * Q2, -j21 * Q2, -(j22 + S) * Q2 - kappa * Q4, dt
    )
    fac = -Q2 * dealias
    p11 *= fac
    p12 *= fac
    p21 *= fac
    p22 *= fac

    m1_0 = 2.0 * (p1b - c11) * (p1b - c12) * (2.0 * p1b - c11 - c12) + (chi + alpha) * p2b + 2.0 * chip * p1b * p2b * p2b
    m2_0 = 2.0 * (p2b - c21) * (p2b - c22) * (2.0 * p2b - c21 - c22) + (chi - alpha) * p1b + 2.0 * chip * p2b * p1b * p1b

    j_acc = np.zeros(len(ts_target), dtype=float)
    for k_seed in range(n_seeds):
        rng = np.random.default_rng(seed0 + k_seed)
        u1 = init_amp * rng.standard_normal((N, N))
        u2 = init_amp * rng.standard_normal((N, N))
        u1 -= u1.mean()
        u2 -= u2.mean()
        u1h = np.fft.rfft2(u1)
        u2h = np.fft.rfft2(u2)

        for step in range(1, nsteps + 1):
            p1 = p1b + u1
            p2 = p2b + u2
            m1 = 2.0 * (p1 - c11) * (p1 - c12) * (2.0 * p1 - c11 - c12) + (chi + alpha) * p2 + 2.0 * chip * p1 * p2 * p2
            m2 = 2.0 * (p2 - c21) * (p2 - c22) * (2.0 * p2 - c21 - c22) + (chi - alpha) * p1 + 2.0 * chip * p2 * p1 * p1
            g1h = np.fft.rfft2((m1 - m1_0) - (j11 + S) * u1 - j12 * u2)
            g2h = np.fft.rfft2((m2 - m2_0) - j21 * u1 - (j22 + S) * u2)
            u1h_new = e11 * u1h + e12 * u2h + p11 * g1h + p12 * g2h
            u2h_new = e21 * u1h + e22 * u2h + p21 * g1h + p22 * g2h
            u1h, u2h = u1h_new, u2h_new
            u1 = np.fft.irfft2(u1h, s=(N, N))
            u2 = np.fft.irfft2(u2h, s=(N, N))

            if step in step_map:
                du1x = np.fft.irfft2(1j * QX * u1h, s=(N, N))
                du1y = np.fft.irfft2(1j * QY * u1h, s=(N, N))
                du2x = np.fft.irfft2(1j * QX * u2h, s=(N, N))
                du2y = np.fft.irfft2(1j * QY * u2h, s=(N, N))
                Jx = float(np.mean(u1 * du2x - u2 * du1x))
                Jy = float(np.mean(u1 * du2y - u2 * du1y))
                j_comp = 0.5 * (abs(Jx) + abs(Jy))
                for idx in step_map[step]:
                    j_acc[idx] += j_comp / n_seeds

    return {
        "J_curve": {"x": ts_target, "y": j_acc.tolist()},
        "J_final": float(j_acc[-1]),
    }
