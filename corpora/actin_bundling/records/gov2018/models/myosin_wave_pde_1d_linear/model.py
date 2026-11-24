import numpy as np


def _subgrid_peak(x_arr, y_arr):
    i = int(np.argmax(y_arr))
    if 0 < i < len(y_arr) - 1:
        y0, y1, y2 = float(y_arr[i - 1]), float(y_arr[i]), float(y_arr[i + 1])
        denom = y0 - 2.0 * y1 + y2
        if abs(denom) > 1e-30:
            delta = 0.5 * (y0 - y2) / denom
            dx = float(x_arr[1] - x_arr[0])
            return float(x_arr[i]) + delta * dx
    return float(x_arr[i])


def run(inputs, settings):
    kappa = float(inputs["kappa"])
    sigma = float(inputs["sigma"])
    kBT = float(inputs["kBT"])
    Lambda = float(inputs["Lambda"])
    H0 = float(inputs["H0"])
    ns = float(inputs["ns"])
    n0 = float(inputs["n0"])
    J = float(inputs["J"])
    xi = float(inputs["xi"])
    A = float(inputs["A"])
    A_star = float(inputs["A_star"])
    k_on = float(inputs["k_on"])
    k_off = float(inputs["k_off"])
    L = float(inputs["L"])
    h0_dip = float(inputs["h0_dip"])
    w0_dip = float(inputs["w0_dip"])
    t_end = float(inputs["t_end"])

    nx = int(round(settings.get("nx", 256)))
    dt = float(settings.get("dt", 0.05))
    n_steps = max(1, int(round(t_end / dt)))
    dt = t_end / n_steps

    x = np.linspace(0.0, L, nx, endpoint=False)
    k = 2.0 * np.pi * np.fft.rfftfreq(nx, d=L / nx)
    k2 = k * k
    k4 = k2 * k2

    h = h0_dip * np.exp(-0.5 * ((x - 0.5 * L) / w0_dip) ** 2)
    h -= np.mean(h)
    dn = np.zeros(nx)
    dm = np.zeros(nx)

    Y = np.stack([np.fft.rfft(h), np.fft.rfft(dn), np.fft.rfft(dm)], axis=1)
    nk = len(k)
    I3 = np.eye(3)
    D = Lambda * kBT
    c_disp = kappa * Lambda * (H0 ** 2) * n0 / ns
    c_curv = -(kappa * Lambda * H0 * n0 / ns)

    A_step = np.zeros((nk, 3, 3), dtype=float)
    for i in range(nk):
        ki2 = k2[i]
        ki4 = k4[i]
        Gh = (kappa * ki4 + sigma * ki2) / xi
        Gn = (D + c_disp) * ki2 - (2.0 * J * Lambda * n0 / ns) * ki2 * (1.0 - ki2 / ns)
        M = np.array(
            [
                [-Gh, (A - kappa * H0 * ki2) / xi, -A_star / xi],
                [c_curv * ki4, -Gn, 0.0],
                [0.0, k_on, -k_off],
            ],
            dtype=float,
        )
        A_step[i] = np.linalg.solve(I3 - 0.5 * dt * M, I3 + 0.5 * dt * M)

    for _ in range(n_steps):
        Y = np.einsum("nij,nj->ni", A_step, Y)

    h_out = np.fft.irfft(Y[:, 0], n=nx)
    n_out = np.fft.irfft(Y[:, 1], n=nx)
    m_out = np.fft.irfft(Y[:, 2], n=nx)

    right_mask = x > 0.5 * L
    x_r = x[right_mask] - 0.5 * L
    x_n = _subgrid_peak(x_r, n_out[right_mask])
    x_m = _subgrid_peak(x_r, m_out[right_mask])
    return {
        "actin_myosin_lead": float(x_n - x_m),
        "wave_packet_speed": float(x_n / t_end),
        "h_min": float(np.min(h_out)),
        "dn_max": float(np.max(n_out)),
        "dn_min": float(np.min(n_out)),
        "dm_max": float(np.max(m_out)),
    }
