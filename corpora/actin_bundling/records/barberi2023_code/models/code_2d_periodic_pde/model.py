import math
import numpy as np
from scipy.fft import rfft2, irfft2
from scipy.optimize import brentq


def run(inputs, settings):
    """Exact line-by-line port of 2D/kernels.jl:3-136 and 2D/main.jl:10-97.
    Uses rfft along axis 0 (x) and fft along axis 1 (y), exact 1-based Gaussian IC
    exp(-((i*dx - L/2)^2 + ...)/(2*Var^2)), traceless velocity denominator
    (1 + visc_coeff*(kx^2 + ky^2)) with visc_coeff = 2 in 2D/kernels.jl:54-55,
    and explicit EulerForward! time stepping (2D/kernels.jl:118-136).
    """
    L = float(inputs["L"])
    Dc = float(inputs["Dc"])
    Da = float(inputs["Da"])
    Di = float(inputs["Di"])
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"])
    Om = float(inputs["Omega"])
    Omd = float(inputs["Omega_d"])
    Omd0 = float(inputs["Omega_d0"])
    Z = float(inputs["Z"])
    B = float(inputs["B"])
    visc_coeff = float(inputs["visc_coeff"])
    var_ic = float(inputs["var_ic"])
    T_end = float(inputs["T_end"])

    N = int(round(float(settings.get("N_grid", 64))))
    dt = float(settings.get("dt", 0.005))
    ic_type = str(settings.get("ic_type", "bump")).lower()
    noise_amp = float(settings.get("noise_amp", 0.0))
    seed = int(round(float(settings.get("seed", 123))))

    dx = L / N
    dy = L / N
    Dap = Da + Di
    Dam = Da - Di

    # 2D/InputParameters.jl:11-14: rfftfreq along x (axis 0), fftfreq along y (axis 1)
    kx = 2.0 * math.pi * np.fft.rfftfreq(N, d=dx)
    ky = 2.0 * math.pi * np.fft.fftfreq(N, d=dy)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K2 = KX * KX + KY * KY
    fac_dx = 1j * KX
    fac_dy = 1j * KY
    fac_lap = -K2
    vel_x_op = (1j * KX) / (1.0 + visc_coeff * K2)
    vel_y_op = (1j * KY) / (1.0 + visc_coeff * K2)

    # 2D/main.jl:49-53: Homogeneous Steady State
    NaH = float(
        brentq(
            lambda Na: -(Omd0 + Omd * A * Na / Kd) * Na + Om0 * (1.0 - Na) * (1.0 + Om * Na * Na),
            0.0,
            1.0,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    CH = (A / Kd) * NaH

    # 2D/kernels.jl:3-29: 1-based indices i = 1..N, j = 1..N
    ii = np.arange(1, N + 1, dtype=float)[:, None]
    jj = np.arange(1, N + 1, dtype=float)[None, :]
    if ic_type == "bump":
        pert = np.exp(-(((ii * dx - 0.5 * L) ** 2 + (jj * dy - 0.5 * L) ** 2) / (2.0 * var_ic * var_ic)))
    elif ic_type == "stripe":
        pert = np.exp(-(((ii * dx - 0.5 * L) ** 2) / (2.0 * var_ic * var_ic))) * np.ones((1, N), dtype=float)
    else:
        pert = np.zeros((N, N), dtype=float)
    pert -= float(np.sum(pert)) * dx * dy / (L * L)

    if noise_amp > 0.0:
        rng = np.random.default_rng(seed)
        noise = noise_amp * (-1.0 + 2.0 * rng.random((N, N)))
        noise -= float(np.sum(noise)) * dx * dy / (L * L)
        pert = pert + noise

    # 2D/main.jl:70-72
    C = CH * (1.0 + pert)
    Nm = (2.0 * NaH - 1.0) * (1.0 + pert)
    Np = 1.0 + pert

    # Guard Euler stability when refining N_grid (keep dt <= 1.0 / (Di * max(K2)))
    max_k2 = float(np.max(K2))
    dt_eff = min(dt, 1.0 / max(1.0, Di * max_k2))
    nsteps = max(1, int(round(T_end / dt_eff)))
    dt_step = T_end / nsteps

    # Exact Fourier-space accumulation of 2D/kernels.jl:43-81,118-136 (linearly identical to real-space Euler)
    fC = rfft2(C, axes=(1, 0))
    fNp = rfft2(Np, axes=(1, 0))
    fNm = rfft2(Nm, axes=(1, 0))
    s_shape = (N, N)

    for _ in range(nsteps):
        Pi = Z * (C * C) - B * (C * C * C)
        fPi = rfft2(Pi, axes=(1, 0))
        Vx = irfft2(fPi * vel_x_op, s=s_shape, axes=(1, 0))
        Vy = irfft2(fPi * vel_y_op, s=s_shape, axes=(1, 0))

        fdiv_CV = fac_dx * rfft2(C * Vx, axes=(1, 0)) + fac_dy * rfft2(C * Vy, axes=(1, 0))
        fdiv_NpV = fac_dx * rfft2(Np * Vx, axes=(1, 0)) + fac_dy * rfft2(Np * Vy, axes=(1, 0))
        fdiv_NmV = fac_dx * rfft2(Nm * Vx, axes=(1, 0)) + fac_dy * rfft2(Nm * Vy, axes=(1, 0))

        fR_Nm = rfft2(
            Om0 * (1.0 + Om * 0.25 * ((Np + Nm) ** 2)) * (Np - Nm) - (Omd0 + Omd * C) * (Np + Nm),
            axes=(1, 0),
        )

        fC_n = fC + dt_step * (-fdiv_CV + (Dc * fac_lap - Kd) * fC + (0.5 * A) * (fNp + fNm))
        fNp_n = fNp + dt_step * (-fdiv_NpV + 0.5 * fac_lap * (Dap * fNp + Dam * fNm))
        fNm_n = fNm + dt_step * (-fdiv_NmV + 0.5 * fac_lap * (Dam * fNp + Dap * fNm) + fR_Nm)
        fC, fNp, fNm = fC_n, fNp_n, fNm_n

        C = irfft2(fC, s=s_shape, axes=(1, 0))
        Np = irfft2(fNp, s=s_shape, axes=(1, 0))
        Nm = irfft2(fNm, s=s_shape, axes=(1, 0))

    return {
        "C_HSS": CH,
        "C_min": float(np.min(C)),
        "C_max": float(np.max(C)),
        "C_corner": float(C[0, 0]),
        "transverse_std": float(np.std(C[N // 2, :])),
        "mean_Np": float(np.mean(Np)),
    }
