import math
import numpy as np


def run(inputs, settings):
    P = float(inputs["P"])
    alpha = float(inputs["alpha"])
    xi = float(inputs["xi"])
    f_0 = float(inputs["f_0"])
    L = float(inputs["L_over_l"])
    t_max = float(inputs["t_max"])
    eps = float(inputs["eps_init"])
    beta_list = [float(b) for b in inputs["T_ref"]["x"]]

    N = int(settings.get("n_grid", 64))
    dt = float(settings.get("dt", 0.03))

    betas = np.array(beta_list, dtype=float)[:, None]
    B = len(beta_list)
    dx = L / N
    x = np.arange(N, dtype=float) * dx
    k = 2.0 * math.pi * np.fft.rfftfreq(N, d=dx)
    k2 = k * k
    dealias = np.abs(k) <= (2.0 / 3.0) * (math.pi / dx)
    scale_f = ((1.0 + xi) ** 2) / xi
    ik = 1j * k
    vel_op = P * ik / (1.0 + k2)
    ik_d = -ik * dealias

    A = np.broadcast_to(1.0 + eps * np.cos(2.0 * math.pi * x / L), (B, N)).copy()
    I = np.broadcast_to(1.0 + eps * np.cos(2.0 * math.pi * x / L), (B, N)).copy()
    Ah = np.fft.rfft(A, axis=1)
    Ih = np.fft.rfft(I, axis=1)

    eA = np.exp(-k2 * dt)
    eI = np.exp(-alpha * k2 * dt)
    phi1A = np.ones_like(k2)
    phi1I = np.ones_like(k2)
    nz = k2 > 0.0
    phi1A[nz] = (1.0 - eA[nz]) / (k2[nz] * dt)
    phi1I[nz] = (1.0 - eI[nz]) / (alpha * k2[nz] * dt)

    n_steps = int(round(t_max / dt))
    half = n_steps // 2
    ts, I_maxs = [], []

    for step in range(n_steps + 1):
        Ar = np.fft.irfft(Ah, n=N, axis=1)
        Ir = np.fft.irfft(Ih, n=N, axis=1)
        fr = scale_f * (f_0 + (1.0 + betas) * Ar / (Ar + xi) + (1.0 - betas) * Ir / (Ir + xi))
        vr = np.fft.irfft(vel_op * np.fft.rfft(fr, axis=1), n=N, axis=1)
        nlA1 = ik_d * np.fft.rfft(vr * Ar, axis=1)
        nlI1 = ik_d * np.fft.rfft(vr * Ir, axis=1)

        if step >= half and step % 2 == 0:
            ts.append(step * dt)
            I_maxs.append(np.max(Ir, axis=1))
        if step == n_steps:
            break

        Ah_p = eA * Ah + dt * phi1A * nlA1
        Ih_p = eI * Ih + dt * phi1I * nlI1
        Ar2 = np.fft.irfft(Ah_p, n=N, axis=1)
        Ir2 = np.fft.irfft(Ih_p, n=N, axis=1)
        fr2 = scale_f * (f_0 + (1.0 + betas) * Ar2 / (Ar2 + xi) + (1.0 - betas) * Ir2 / (Ir2 + xi))
        vr2 = np.fft.irfft(vel_op * np.fft.rfft(fr2, axis=1), n=N, axis=1)
        nlA2 = ik_d * np.fft.rfft(vr2 * Ar2, axis=1)
        nlI2 = ik_d * np.fft.rfft(vr2 * Ir2, axis=1)

        Ah = eA * Ah + 0.5 * dt * (eA * nlA1 + nlA2)
        Ih = eI * Ih + 0.5 * dt * (eI * nlI1 + nlI2)

    ts = np.array(ts)
    I_maxs = np.array(I_maxs)
    Ts = []
    for j in range(B):
        Ih_ = I_maxs[:, j]
        osc_amp = float(np.max(Ih_) - np.min(Ih_))
        if osc_amp < 5e-2:
            Ts.append(0.0)
            continue
        pks = np.where((Ih_[1:-1] > Ih_[:-2]) & (Ih_[1:-1] >= Ih_[2:]) & (Ih_[1:-1] > 1.15))[0] + 1
        t_pks = []
        for p in pks:
            y0, y1, y2 = Ih_[p - 1], Ih_[p], Ih_[p + 1]
            denom = y0 - 2.0 * y1 + y2
            delta = 0.5 * (y0 - y2) / denom if denom != 0.0 else 0.0
            t_pks.append(ts[p] + delta * (ts[p + 1] - ts[p]))
        if len(t_pks) >= 3:
            Ts.append(float(2.0 * np.mean(np.diff(t_pks))))
        elif len(t_pks) == 2:
            Ts.append(float(2.0 * (t_pks[1] - t_pks[0])))
        else:
            Ts.append(0.0)
    return {"T_over_tau": {"x": beta_list, "y": Ts}}
