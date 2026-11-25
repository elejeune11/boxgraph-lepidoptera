"""2D pseudospectral IMEX Cahn-Hilliard solver returning mean order parameter <phi>(t) on the Fig. S6 time grid."""

import numpy as np
from scipy import fft as spfft


def run(inputs, settings):
    M = float(inputs["M"])
    kappa = float(inputs["kappa"])
    gamma_dot = float(inputs["gamma_dot"])
    L = float(inputs["L"])
    phi_avg = float(inputs["phi_avg"])
    noise_amp = float(inputs["noise_amp"])
    t_targets = [float(t) for t in inputs["ref_mass"]["x"]]

    flow_type = str(settings.get("flow_type", "rotational"))
    N = int(settings.get("N", 128))
    dt = float(settings.get("dt", 0.05))
    seed = int(settings.get("seed", 0))

    rng = np.random.default_rng(seed)
    dx = L / N
    coords = (np.arange(N) - 0.5 * N) * dx
    X, Y = np.meshgrid(coords, coords, indexing="xy")
    # Zero-mean noise preserves exact initial <phi>(0) = phi_avg
    noise = rng.standard_normal((N, N))
    noise -= float(np.mean(noise))
    phi = phi_avg + noise_amp * noise

    kx = spfft.rfftfreq(N, d=dx) * 2.0 * np.pi
    ky = spfft.fftfreq(N, d=dx) * 2.0 * np.pi
    KX, KY = np.meshgrid(kx, ky, indexing="xy")
    K2 = KX * KX + KY * KY
    K4 = K2 * K2
    k_nyq = np.pi / dx
    dealias = (np.abs(KX) <= (2.0 / 3.0) * k_nyq) & (np.abs(KY) <= (2.0 / 3.0) * k_nyq)

    if flow_type == "shear":
        vx = gamma_dot * Y
        vy = np.zeros_like(X)
    elif flow_type == "rotational":
        vx = gamma_dot * Y
        vy = -gamma_dot * X
    else:
        raise ValueError(f"unknown flow_type: {flow_type!r}")

    phi_k = spfft.rfft2(phi)
    denom = 1.0 + dt * M * kappa * K4
    target_steps = {int(round(t / dt)): idx for idx, t in enumerate(t_targets)}
    nsteps = max(target_steps.keys())
    out = np.zeros(len(t_targets), dtype=np.float64)

    for step in range(1, nsteps + 1):
        phi_cur = spfft.irfft2(phi_k, s=(N, N))
        df_k = spfft.rfft2(phi_cur * phi_cur * phi_cur - phi_cur)
        df_k *= dealias
        dphidx = spfft.irfft2(1j * KX * phi_k, s=(N, N))
        dphidy = spfft.irfft2(1j * KY * phi_k, s=(N, N))
        adv_k = spfft.rfft2(vx * dphidx + vy * dphidy)
        adv_k *= dealias
        phi_k = (phi_k + dt * (-M * K2 * df_k - adv_k)) / denom
        if step in target_steps:
            out[target_steps[step]] = float(np.mean(spfft.irfft2(phi_k, s=(N, N))))

    return {"mass": {"x": t_targets, "y": [float(v) for v in out]}, "mean_phi_final": float(out[-1])}
