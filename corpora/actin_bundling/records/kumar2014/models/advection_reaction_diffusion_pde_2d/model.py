import numpy as np


def run(inputs, settings):
    P = float(inputs["P"])
    beta = float(inputs["beta"])
    K = float(inputs["kappa_tau"])
    rho = float(inputs["rho"])
    nu = float(inputs["nu"])
    xi = float(inputs["xi"])
    f_0 = float(inputs["f_0"])
    L = float(inputs["L_over_l"])
    t_max = float(inputs["t_max"])
    eps = float(inputs["eps_init"])

    N = int(settings.get("n_grid", 64))
    dt = float(settings.get("dt", 0.04))
    seed = int(settings.get("seed", 42))

    rng = np.random.default_rng(seed)
    dx = L / N
    coords = np.arange(N, dtype=float) * dx
    X, Y = np.meshgrid(coords, coords, indexing="ij")

    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    ky = 2.0 * np.pi * np.fft.rfftfreq(N, d=dx)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K2 = KX * KX + KY * KY
    kmax = (2.0 / 3.0) * (np.pi / dx)
    dealias = (np.abs(KX) <= kmax) & (np.abs(KY) <= kmax)
    scale_f = ((1.0 + xi) ** 2) / xi

    # Fixed low-mode seeded 2D Fourier perturbation so grid refinement preserves the exact IC
    A = np.ones((N, N), dtype=float)
    I = np.ones((N, N), dtype=float)
    modes = [(mx, my) for mx in range(-4, 5) for my in range(-4, 5) if 1 <= mx * mx + my * my <= 16]
    norm = eps / np.sqrt(len(modes))
    for mx, my in modes:
        phaseA = rng.uniform(0.0, 2.0 * np.pi)
        phaseI = rng.uniform(0.0, 2.0 * np.pi)
        ampA = norm * rng.uniform(0.5, 1.0)
        ampI = norm * rng.uniform(0.5, 1.0)
        arg = 2.0 * np.pi * (mx * X + my * Y) / L
        A += ampA * np.cos(arg + phaseA)
        I += ampI * np.cos(arg + phaseI)

    Ah = np.fft.rfft2(A)
    Ih = np.fft.rfft2(I)

    lamA = K2 + K
    lamI = K2 + rho * K
    eA = np.exp(-lamA * dt)
    eI = np.exp(-lamI * dt)
    phiA = (1.0 - eA) / (lamA * dt)
    phiI = (1.0 - eI) / (lamI * dt)

    srcA = np.zeros_like(Ah)
    srcI = np.zeros_like(Ih)
    srcA[0, 0] = K * (N * N)
    srcI[0, 0] = rho * K * (N * N)

    denom_v = 1.0 + nu * K2

    def rhs_nl(Ah_, Ih_):
        Ar = np.fft.irfft2(Ah_, s=(N, N))
        Ir = np.fft.irfft2(Ih_, s=(N, N))
        fr = scale_f * (f_0 + (1.0 + beta) * Ar / (Ar + xi) + (1.0 - beta) * Ir / (Ir + xi))
        fh = np.fft.rfft2(fr)
        vx = np.fft.irfft2(P * 1j * KX * fh / denom_v, s=(N, N))
        vy = np.fft.irfft2(P * 1j * KY * fh / denom_v, s=(N, N))
        nlA = (-1j * KX * np.fft.rfft2(vx * Ar) - 1j * KY * np.fft.rfft2(vy * Ar)) * dealias + srcA
        nlI = (-1j * KX * np.fft.rfft2(vx * Ir) - 1j * KY * np.fft.rfft2(vy * Ir)) * dealias + srcI
        return nlA, nlI, Ar, Ir

    n_steps = int(round(t_max / dt))
    sample_every = max(1, int(round(0.5 / dt)))
    A_maxs, I_maxs, A_pt0 = [], [], []

    for step in range(n_steps + 1):
        nlA1, nlI1, Ar, Ir = rhs_nl(Ah, Ih)
        if step % sample_every == 0:
            A_maxs.append(float(np.max(Ar)))
            I_maxs.append(float(np.max(Ir)))
            A_pt0.append(float(Ar[0, 0]))
        if step == n_steps:
            break
        Ah_p = eA * Ah + dt * phiA * nlA1
        Ih_p = eI * Ih + dt * phiI * nlI1
        nlA2, nlI2, _, _ = rhs_nl(Ah_p, Ih_p)
        Ah = eA * Ah + 0.5 * dt * (eA * nlA1 + nlA2)
        Ih = eI * Ih + 0.5 * dt * (eI * nlI1 + nlI2)

    Ar_final = np.fft.irfft2(Ah, s=(N, N))
    Ir_final = np.fft.irfft2(Ih, s=(N, N))
    A_maxs = np.array(A_maxs)
    I_maxs = np.array(I_maxs)
    A_pt0 = np.array(A_pt0)

    # Evaluate post-transient last third
    start = (2 * len(A_maxs)) // 3
    A_tail = A_maxs[start:]
    I_tail = I_maxs[start:]
    pt_tail = A_pt0[start:]

    spatial_dev = float(np.max(A_tail) - 1.0)
    osc_amp = float(np.max(A_tail) - np.min(A_tail))
    pt_amp = float(np.max(pt_tail) - np.min(pt_tail))

    if spatial_dev < 2e-2:
        state_code = 0.0
    elif osc_amp > 0.1 or pt_amp > 0.1:
        state_code = 2.0
    else:
        state_code = 1.0

    # Count distinct local maxima above 1.2 in final A field
    is_max = np.ones((N, N), dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            is_max &= Ar_final >= np.roll(np.roll(Ar_final, di, axis=0), dj, axis=1)
    n_peaks = float(np.sum(is_max & (Ar_final > 1.2)))

    return {
        "state_code": state_code,
        "A_max": float(np.max(A_tail)),
        "I_max": float(np.max(I_tail)),
        "osc_amp": osc_amp,
        "n_peaks": n_peaks,
        "A_mean_final": float(np.mean(Ar_final)),
        "I_mean_final": float(np.mean(Ir_final)),
    }
