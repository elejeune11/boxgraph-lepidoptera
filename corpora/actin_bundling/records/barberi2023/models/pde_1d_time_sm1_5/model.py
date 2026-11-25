import math
import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    L = float(inputs["L"])
    Dc = float(inputs["Dc"])
    Da = float(inputs["Da"])
    Di = float(inputs["Di"])
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"])
    Om = float(inputs["Omega"])
    Omd0 = float(inputs["Omega_d0"])
    Omd = float(inputs["Omega_d"])
    Z = float(inputs["Z"])
    B = float(inputs["B"])
    Nbar = float(inputs["Nbar"])
    T_end = float(inputs["T_end"])

    N = int(round(settings.get("n_grid", 128)))
    dt = float(settings.get("dt", 0.05))
    noise_amp = float(settings.get("noise_amp", 0.01))
    seed = int(round(settings.get("seed", 42)))

    ratio = A / Kd
    na_h = float(
        brentq(
            lambda x: Om0 * (1.0 + Om * x * x) * (Nbar - x) - (Omd0 + Omd * ratio * x) * x,
            0.0,
            Nbar,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    c_h = ratio * na_h
    ni_h = Nbar - na_h

    rng = np.random.default_rng(seed)
    C = c_h + noise_amp * rng.standard_normal(N)
    Na = na_h + noise_amp * rng.standard_normal(N)
    Ni = ni_h + noise_amp * rng.standard_normal(N)
    Ni += Nbar - float(np.mean(Na + Ni))

    k = 2.0 * math.pi * np.fft.rfftfreq(N, d=L / N)
    k2 = k * k
    ik = 1j * k
    v_op = ik / (1.0 + k2)
    cutoff = (len(k) * 2) // 3

    nsteps = max(1, int(round(T_end / dt)))
    dt_act = T_end / nsteps
    eC = np.exp(-(Dc * k2 + Kd) * dt_act)
    eNa = np.exp(-Da * k2 * dt_act)
    eNi = np.exp(-Di * k2 * dt_act)

    Ch = np.fft.rfft(C)
    Nah = np.fft.rfft(Na)
    Nih = np.fft.rfft(Ni)
    active = (Z != 0.0) or (B != 0.0)

    def rhs(c_h_arr, na_h_arr, ni_h_arr):
        c_cur = np.fft.irfft(c_h_arr, n=N)
        na_cur = np.fft.irfft(na_h_arr, n=N)
        ni_cur = np.fft.irfft(ni_h_arr, n=N)
        r_h = np.fft.rfft(Om0 * (1.0 + Om * na_cur * na_cur) * ni_cur - (Omd0 + Omd * c_cur) * na_cur)
        if active:
            pi_h = np.fft.rfft(Z * c_cur * c_cur - B * c_cur * c_cur * c_cur)
            v_cur = np.fft.irfft(v_op * pi_h, n=N)
            nc = -ik * np.fft.rfft(v_cur * c_cur) + A * na_h_arr
            nna = -ik * np.fft.rfft(v_cur * na_cur) + r_h
            nni = -ik * np.fft.rfft(v_cur * ni_cur) - r_h
        else:
            nc = A * na_h_arr
            nna = r_h
            nni = -r_h
        nc[cutoff:] = 0.0
        nna[cutoff:] = 0.0
        nni[cutoff:] = 0.0
        return nc, nna, nni

    for _ in range(nsteps):
        nc1, nna1, nni1 = rhs(Ch, Nah, Nih)
        cp = eC * (Ch + dt_act * nc1)
        nap = eNa * (Nah + dt_act * nna1)
        nip = eNi * (Nih + dt_act * nni1)
        nc2, nna2, nni2 = rhs(cp, nap, nip)
        Ch = eC * Ch + 0.5 * dt_act * (eC * nc1 + nc2)
        Nah = eNa * Nah + 0.5 * dt_act * (eNa * nna1 + nna2)
        Nih = eNi * Nih + 0.5 * dt_act * (eNi * nni1 + nni2)

    C_f = np.fft.irfft(Ch, n=N)
    amp = np.abs(np.fft.rfft(C_f))[1:]
    n_dom = int(np.argmax(amp)) + 1
    k_dom = 2.0 * math.pi * n_dom / L

    return {
        "C_min": float(np.min(C_f)),
        "C_max": float(np.max(C_f)),
        "K_dom": float(k_dom),
    }
