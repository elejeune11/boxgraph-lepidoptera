import math
import numpy as np


def simulate_1d(inputs, settings):
    P = float(inputs["P"])
    beta = float(inputs["beta"])
    alpha = float(inputs["alpha"])
    xi = float(inputs["xi"])
    f_0 = float(inputs["f_0"])
    L = float(inputs["L_over_l"])
    t_max = float(inputs["t_max"])
    eps = float(inputs["eps_init"])

    N = int(settings.get("n_grid", 64))
    dt = float(settings.get("dt", 0.01))
    seed = int(settings.get("seed", 42))

    rng = np.random.default_rng(seed)
    dx = L / N
    x = np.arange(N, dtype=float) * dx
    k = 2.0 * np.pi * np.fft.rfftfreq(N, d=dx)
    k2 = k * k
    dealias = np.abs(k) <= (2.0 / 3.0) * (np.pi / dx)
    scale_f = ((1.0 + xi) ** 2) / xi

    # Fixed low-mode seeded random perturbation so grid refinement preserves the exact IC
    A = np.ones(N, dtype=float)
    I = np.ones(N, dtype=float)
    for m in range(1, 5):
        ampA = eps * (0.5 ** (m - 1)) * rng.uniform(0.5, 1.0)
        ampI = eps * (0.5 ** (m - 1)) * rng.uniform(0.5, 1.0)
        phiA = rng.uniform(0.0, 2.0 * np.pi)
        phiI = rng.uniform(0.0, 2.0 * np.pi)
        A += ampA * np.cos(2.0 * np.pi * m * x / L + phiA)
        I += ampI * np.cos(2.0 * np.pi * m * x / L + phiI)

    Ah = np.fft.rfft(A)
    Ih = np.fft.rfft(I)

    eA = np.exp(-k2 * dt)
    eI = np.exp(-alpha * k2 * dt)
    phi1A = np.ones_like(k2)
    phi1I = np.ones_like(k2)
    nz = k2 > 0.0
    phi1A[nz] = (1.0 - eA[nz]) / (k2[nz] * dt)
    phi1I[nz] = (1.0 - eI[nz]) / (alpha * k2[nz] * dt)

    def rhs_nl(Ah_, Ih_):
        Ar = np.fft.irfft(Ah_, n=N)
        Ir = np.fft.irfft(Ih_, n=N)
        fr = scale_f * (f_0 + (1.0 + beta) * Ar / (Ar + xi) + (1.0 - beta) * Ir / (Ir + xi))
        fh = np.fft.rfft(fr)
        vh = P * 1j * k * fh / (1.0 + k2)
        vr = np.fft.irfft(vh, n=N)
        nlA = -1j * k * np.fft.rfft(vr * Ar) * dealias
        nlI = -1j * k * np.fft.rfft(vr * Ir) * dealias
        return nlA, nlI, Ar, Ir, vr, fr

    n_steps = int(round(t_max / dt))
    sample_dt = 0.04
    sample_every = max(1, int(round(sample_dt / dt)))

    ts, A_maxs, I_maxs, v_maxs, f_maxs, f_mins = [], [], [], [], [], []
    A_pt0 = []
    for step in range(n_steps + 1):
        nlA1, nlI1, Ar, Ir, vr, fr = rhs_nl(Ah, Ih)
        if step % sample_every == 0:
            ts.append(step * dt)
            A_maxs.append(float(np.max(Ar)))
            I_maxs.append(float(np.max(Ir)))
            v_maxs.append(float(np.max(np.abs(vr))))
            f_maxs.append(float(np.max(fr)))
            f_mins.append(float(np.min(fr)))
            A_pt0.append(float(Ar[0]))
        if step == n_steps:
            break
        Ah_p = eA * Ah + dt * phi1A * nlA1
        Ih_p = eI * Ih + dt * phi1I * nlI1
        nlA2, nlI2, _, _, _, _ = rhs_nl(Ah_p, Ih_p)
        Ah = eA * Ah + 0.5 * dt * (eA * nlA1 + nlA2)
        Ih = eI * Ih + 0.5 * dt * (eI * nlI1 + nlI2)

    Ar_final = np.fft.irfft(Ah, n=N)
    Ir_final = np.fft.irfft(Ih, n=N)
    ts = np.array(ts)
    A_maxs = np.array(A_maxs)
    I_maxs = np.array(I_maxs)
    v_maxs = np.array(v_maxs)
    f_maxs = np.array(f_maxs)
    f_mins = np.array(f_mins)
    A_pt0 = np.array(A_pt0)

    half = len(ts) // 2
    I_h = I_maxs[half:]
    A_h = A_maxs[half:]
    t_h = ts[half:]

    spatial_dev = max(float(np.max(I_h)) - 1.0, float(np.max(A_h)) - 1.0)
    osc_amp = float(np.max(I_h) - np.min(I_h))
    pt_amp = float(np.max(A_pt0[half:]) - np.min(A_pt0[half:]))

    if spatial_dev < 1e-2:
        state_code = 0.0
        T_osc = 0.0
    elif osc_amp > 5e-2 or pt_amp > 5e-2:
        state_code = 2.0
        # Standing-wave pulse flips by L/2 every T/2 -> peaks of max(I) occur every T/2
        pks = np.where((I_h[1:-1] > I_h[:-2]) & (I_h[1:-1] >= I_h[2:]) & (I_h[1:-1] > 1.15))[0] + 1
        t_pks = []
        for p in pks:
            y0, y1, y2 = I_h[p - 1], I_h[p], I_h[p + 1]
            denom = y0 - 2.0 * y1 + y2
            delta = 0.5 * (y0 - y2) / denom if denom != 0.0 else 0.0
            t_pks.append(t_h[p] + delta * (t_h[p + 1] - t_h[p]))
        if len(t_pks) >= 3:
            T_osc = float(2.0 * np.mean(np.diff(t_pks)))
        else:
            T_osc = 0.0
    else:
        state_code = 1.0
        T_osc = 0.0

    return {
        "state_code": state_code,
        "T_over_tau": T_osc,
        "A_max": float(np.max(A_h)),
        "I_max": float(np.max(I_h)),
        "v_max": float(np.max(v_maxs[half:])),
        "f_max": float(np.max(f_maxs[half:])),
        "f_min": float(np.min(f_mins[half:])),
        "A_mean_final": float(np.mean(Ar_final)),
        "I_mean_final": float(np.mean(Ir_final)),
        "scale_f": scale_f,
    }


def run(inputs, settings):
    res = simulate_1d(inputs, settings)
    res.pop("scale_f")
    return res
