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

    N = int(round(settings.get("n_grid", 64)))
    dt = float(settings.get("dt", 0.025))
    seed_w = float(settings.get("seed_width", 2.0))

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

    x = np.linspace(-0.5 * L, 0.5 * L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    R2 = X * X + Y * Y
    bump = np.exp(-R2 / (seed_w * seed_w))
    C = c_h + 0.7 * bump
    Na = na_h + 0.6 * bump
    Ni = ni_h + 0.5 * bump
    scale = Nbar / float(np.mean(Na + Ni))
    Na *= scale
    Ni *= scale

    kx = 2.0 * math.pi * np.fft.fftfreq(N, d=L / N)
    ky = 2.0 * math.pi * np.fft.rfftfreq(N, d=L / N)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K2 = KX * KX + KY * KY
    kmax = float(np.max(np.abs(kx)))
    dealias = (np.abs(KX) <= (2.0 / 3.0) * kmax) & (np.abs(KY) <= (2.0 / 3.0) * kmax)

    # Printed Eq. (3) on p. 1 defines v = (1/2)[grad v + (grad v)^T - (div v / d) 1].
    # In d = 2 with eta = eta_tilde: div sigma_visc = eta lap(v) + (3/2) eta grad(div v),
    # which for longitudinal flow gives (1 + 2.5 |K|^2) V_hat = i K Pi_hat.
    vx_op = 1j * KX / (1.0 + 2.5 * K2)
    vy_op = 1j * KY / (1.0 + 2.5 * K2)

    nsteps = max(1, int(round(T_end / dt)))
    dt_act = T_end / nsteps
    inv_C = 1.0 / (1.0 + dt_act * (Dc * K2 + Kd))
    inv_Na = 1.0 / (1.0 + dt_act * Da * K2)
    inv_Ni = 1.0 / (1.0 + dt_act * Di * K2)

    Ch = np.fft.rfft2(C)
    Nah = np.fft.rfft2(Na)
    Nih = np.fft.rfft2(Ni)
    for _ in range(nsteps):
        C = np.fft.irfft2(Ch, s=(N, N))
        Na = np.fft.irfft2(Nah, s=(N, N))
        Ni = np.fft.irfft2(Nih, s=(N, N))
        pih = np.fft.rfft2(Z * C * C - B * C * C * C)
        vx = np.fft.irfft2(vx_op * pih, s=(N, N))
        vy = np.fft.irfft2(vy_op * pih, s=(N, N))
        rh = np.fft.rfft2(Om0 * (1.0 + Om * Na * Na) * Ni - (Omd0 + Omd * C) * Na) * dealias
        advC = -(1j * KX * np.fft.rfft2(vx * C) + 1j * KY * np.fft.rfft2(vy * C)) * dealias
        advNa = -(1j * KX * np.fft.rfft2(vx * Na) + 1j * KY * np.fft.rfft2(vy * Na)) * dealias
        advNi = -(1j * KX * np.fft.rfft2(vx * Ni) + 1j * KY * np.fft.rfft2(vy * Ni)) * dealias
        Ch = (Ch + dt_act * (advC + A * Nah)) * inv_C
        Nah = (Nah + dt_act * (advNa + rh)) * inv_Na
        Nih = (Nih + dt_act * (advNi - rh)) * inv_Ni

    C_f = np.fft.irfft2(Ch, s=(N, N))
    return {
        "C_min": float(np.min(C_f)),
        "C_max": float(np.max(C_f)),
        "C_corner": float(C_f[0, 0]),
    }
