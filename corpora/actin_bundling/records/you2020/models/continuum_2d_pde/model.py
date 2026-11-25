import math
import numpy as np


def run(inputs, settings):
    chiA = float(inputs["chi_A"])
    chiB = float(inputs["chi_B"])
    gammaA = float(inputs["gamma_A"])
    kappa = float(inputs["kappa"])
    delta = float(inputs["delta"])
    L = float(inputs["L"])
    phi0A = float(inputs["phi0_A"])
    phi0B = float(inputs["phi0_B"])

    N = int(round(float(settings.get("N", 32))))
    dt = float(settings.get("dt", 2.0))
    t_max = float(settings.get("t_max", 10000.0))
    init_mode = str(settings.get("init_mode", "stripe"))
    seed = int(settings.get("seed", 7))

    kx = 2.0 * math.pi * np.fft.fftfreq(N, d=L / N)
    ky = 2.0 * math.pi * np.fft.rfftfreq(N, d=L / N)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K2 = KX * KX + KY * KY
    kAB = kappa - delta
    kBA = kappa + delta

    ix = np.abs(np.fft.fftfreq(N) * N)
    iy = np.abs(np.fft.rfftfreq(N) * N)
    IX, IY = np.meshgrid(ix, iy, indexing="ij")
    dealias = (IX <= N // 3) & (IY <= N // 3)

    x = np.linspace(0.0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    rng = np.random.default_rng(seed)

    if init_mode == "oscillatory":
        phiA = phi0A + 0.05 * (np.cos(2.0 * math.pi * X / L) + np.cos(2.0 * math.pi * Y / L)) + 1.0e-4 * rng.standard_normal((N, N))
        phiB = phi0B - 0.05 * (np.cos(2.0 * math.pi * X / L) - np.cos(2.0 * math.pi * Y / L)) + 1.0e-4 * rng.standard_normal((N, N))
    else:
        phiA = phi0A + 1.0e-2 * np.cos(2.0 * math.pi * X / L) + 1.0e-4 * rng.standard_normal((N, N))
        phiB = phi0B - 1.0e-2 * np.cos(2.0 * math.pi * X / L + 0.05) + 1.0e-4 * rng.standard_normal((N, N))

    phiA -= phiA.mean() - phi0A
    phiB -= phiB.mean() - phi0B

    Ahat = np.fft.rfft2(phiA)
    Bhat = np.fft.rfft2(phiB)

    M11 = -K2 * (chiA + gammaA * K2)
    M12 = -K2 * kAB
    M21 = -K2 * kBA
    M22 = -K2 * chiB

    e11 = 1.0 - dt * M11
    e12 = -dt * M12
    e21 = -dt * M21
    e22 = 1.0 - dt * M22
    edet = e11 * e22 - e12 * e21

    s11 = 3.0 - 2.0 * dt * M11
    s12 = -2.0 * dt * M12
    s21 = -2.0 * dt * M21
    s22 = 3.0 - 2.0 * dt * M22
    sdet = s11 * s22 - s12 * s21

    steps = int(round(t_max / dt))
    window_time = min(6000.0, 0.6 * t_max)
    sample_start = max(0, steps - int(round(window_time / dt)))
    sample_every = max(1, int(round(10.0 / dt)))

    A = np.fft.irfft2(Ahat, s=(N, N))
    nlA_prev = -(1.0 / 3.0) * K2 * np.fft.rfft2(A * A * A) * dealias
    Ahat_prev = Ahat.copy()
    Bhat_prev = Bhat.copy()

    rhsA = Ahat + dt * nlA_prev
    rhsB = Bhat
    Ahat = (e22 * rhsA - e12 * rhsB) / edet
    Bhat = (-e21 * rhsA + e11 * rhsB) / edet

    phases = []
    mode_re = []
    times = []
    for n in range(1, steps):
        A = np.fft.irfft2(Ahat, s=(N, N))
        nlA = -(1.0 / 3.0) * K2 * np.fft.rfft2(A * A * A) * dealias

        rhsA = 4.0 * Ahat - Ahat_prev + 2.0 * dt * (2.0 * nlA - nlA_prev)
        rhsB = 4.0 * Bhat - Bhat_prev

        Ahat_prev = Ahat
        Bhat_prev = Bhat
        nlA_prev = nlA

        Ahat = (s22 * rhsA - s12 * rhsB) / sdet
        Bhat = (-s21 * rhsA + s11 * rhsB) / sdet

        if n >= sample_start and (n - sample_start) % sample_every == 0:
            phases.append(math.atan2(Ahat[1, 0].imag, Ahat[1, 0].real))
            mode_re.append(float(Ahat[1, 0].real / (N * N)))
            times.append(n * dt)

    times_arr = np.asarray(times, dtype=float)
    amp1 = float(abs(Ahat[1, 0]) / (N * N))
    amp2 = float(abs(Ahat[0, 1]) / (N * N))
    if max(amp1, amp2) < 1.0e-4:
        return {"v": 0.0, "freq_osc": 0.0, "state_code": 0.0}

    if init_mode == "oscillatory":
        sig = np.asarray(mode_re, dtype=float)
        sig = sig - sig.mean()
        if np.max(np.abs(sig)) < 1.0e-4:
            return {"v": 0.0, "freq_osc": 0.0, "state_code": 1.0}
        up = []
        for i in range(len(sig) - 1):
            if sig[i] <= 0.0 and sig[i + 1] > 0.0:
                frac = -sig[i] / (sig[i + 1] - sig[i])
                up.append(times_arr[i] + frac * (times_arr[i + 1] - times_arr[i]))
        if len(up) >= 2:
            period = (up[-1] - up[0]) / (len(up) - 1)
            freq = float(1.0 / period)
        else:
            spec = np.abs(np.fft.rfft(sig))
            freqs = np.fft.rfftfreq(len(sig), d=times_arr[1] - times_arr[0])
            freq = float(freqs[int(np.argmax(spec[1:])) + 1])
        return {"v": 0.0, "freq_osc": freq, "state_code": 2.0 if freq > 1e-5 else 1.0}

    ph = np.unwrap(phases)
    q1 = 2.0 * math.pi / L
    v = float(abs((ph[-1] - ph[0]) / (times_arr[-1] - times_arr[0])) / q1)
    state = 2.0 if v > 1.0e-4 else 1.0
    return {"v": v, "freq_osc": 0.0, "state_code": state}
