import math
import numpy as np


def run(inputs, settings):
    alpha = float(inputs["alpha"])
    xi = float(inputs["xi"])
    f_0 = float(inputs["f_0"])
    L = float(inputs["L_over_l"])
    t_max = float(inputs["t_max"])
    eps = float(inputs["eps_init"])
    beta_list = [float(b) for b in inputs["bnd_ref"]["x"]]

    N = int(settings.get("n_grid", 48))
    dt = float(settings.get("dt", 0.03))
    n_bisect = int(settings.get("n_bisect", 6))

    betas = np.array(beta_list, dtype=float)[:, None]
    B = len(beta_list)
    dx = L / N
    x = np.arange(N, dtype=float) * dx
    k = 2.0 * math.pi * np.fft.rfftfreq(N, d=dx)
    k2 = k * k
    dealias = np.abs(k) <= (2.0 / 3.0) * (math.pi / dx)
    scale_f = ((1.0 + xi) ** 2) / xi
    ik = 1j * k
    ik_d = -ik * dealias

    eA = np.exp(-k2 * dt)
    eI = np.exp(-alpha * k2 * dt)
    phi1A = np.ones_like(k2)
    phi1I = np.ones_like(k2)
    nz = k2 > 0.0
    phi1A[nz] = (1.0 - eA[nz]) / (k2[nz] * dt)
    phi1I[nz] = (1.0 - eI[nz]) / (alpha * k2[nz] * dt)

    n_steps = int(round(t_max / dt))
    half = n_steps // 2
    A0 = np.broadcast_to(1.0 + eps * np.cos(2.0 * math.pi * x / L), (B, N))

    lo = np.full(B, 1.15, dtype=float)
    hi = np.full(B, 2.05, dtype=float)

    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        vel_op = mid[:, None] * ik / (1.0 + k2)
        Ah = np.fft.rfft(A0.copy(), axis=1)
        Ih = np.fft.rfft(A0.copy(), axis=1)
        I_hi = np.full(B, -1e9)
        I_lo = np.full(B, 1e9)

        for step in range(n_steps + 1):
            Ar = np.fft.irfft(Ah, n=N, axis=1)
            Ir = np.fft.irfft(Ih, n=N, axis=1)
            fr = scale_f * (f_0 + (1.0 + betas) * Ar / (Ar + xi) + (1.0 - betas) * Ir / (Ir + xi))
            vr = np.fft.irfft(vel_op * np.fft.rfft(fr, axis=1), n=N, axis=1)
            nlA1 = ik_d * np.fft.rfft(vr * Ar, axis=1)
            nlI1 = ik_d * np.fft.rfft(vr * Ir, axis=1)

            if step >= half and step % 4 == 0:
                im = np.max(Ir, axis=1)
                I_hi = np.maximum(I_hi, im)
                I_lo = np.minimum(I_lo, im)
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

        is_osc = (I_hi - I_lo) > 0.05
        lo = np.where(is_osc, mid, lo)
        hi = np.where(is_osc, hi, mid)

    P_bnd = 0.5 * (lo + hi)
    return {"P_osc_to_stat": {"x": beta_list, "y": [float(v) for v in P_bnd]}}
