"""2D pseudospectral IMEX Cahn-Hilliard solver returning structure-factor domain size a(t) on the reference t-grid."""

import numpy as np
from scipy import fft as spfft


def _single_at(M, kappa, gamma_dot, L, phi_avg, noise_amp, flow_type, N, dt, seed, t_targets):
    rng = np.random.default_rng(seed)
    dx = L / N
    coords = (np.arange(N) - 0.5 * N) * dx
    X, Y = np.meshgrid(coords, coords, indexing="xy")
    phi = phi_avg + noise_amp * rng.standard_normal((N, N))

    kx = spfft.rfftfreq(N, d=dx) * 2.0 * np.pi
    ky = spfft.fftfreq(N, d=dx) * 2.0 * np.pi
    KX, KY = np.meshgrid(kx, ky, indexing="xy")
    K2 = KX * KX + KY * KY
    K4 = K2 * K2
    k_nyq = np.pi / dx
    dealias = (np.abs(KX) <= (2.0 / 3.0) * k_nyq) & (np.abs(KY) <= (2.0 / 3.0) * k_nyq)

    if flow_type == "passive" or gamma_dot == 0.0:
        vx = vy = None
    elif flow_type == "shear":
        vx = gamma_dot * Y
        vy = np.zeros_like(X)
    elif flow_type == "rotational":
        vx = gamma_dot * Y
        vy = -gamma_dot * X
    else:
        raise ValueError(f"unknown flow_type: {flow_type!r}")

    dk = 2.0 * np.pi / L
    K_mag = np.sqrt(K2)
    bins = np.round(K_mag / dk).astype(np.int64)
    max_b = N // 3
    w_rfft = np.ones_like(KX)
    w_rfft[:, 1:-1] = 2.0
    mask = (bins >= 1) & (bins <= max_b)
    bins_m = bins[mask]
    w_m = w_rfft[mask]
    counts = np.bincount(bins_m, weights=w_m, minlength=max_b + 1)
    ks = np.arange(1, max_b + 1) * dk

    phi_k = spfft.rfft2(phi)
    denom = 1.0 + dt * M * kappa * K4
    target_steps = {int(round(t / dt)): idx for idx, t in enumerate(t_targets)}
    nsteps = max(target_steps.keys())
    out = np.zeros(len(t_targets), dtype=np.float64)

    for step in range(1, nsteps + 1):
        phi_cur = spfft.irfft2(phi_k, s=(N, N))
        df_k = spfft.rfft2(phi_cur * phi_cur * phi_cur - phi_cur)
        df_k *= dealias
        rhs = -M * K2 * df_k
        if vx is not None:
            dphidx = spfft.irfft2(1j * KX * phi_k, s=(N, N))
            dphidy = spfft.irfft2(1j * KY * phi_k, s=(N, N))
            adv_k = spfft.rfft2(vx * dphidx + vy * dphidy)
            adv_k *= dealias
            rhs -= adv_k
        phi_k = (phi_k + dt * rhs) / denom

        if step in target_steps:
            S2d = np.abs(phi_k) ** 2
            S_sum = np.bincount(bins_m, weights=(S2d[mask] * w_m), minlength=max_b + 1)
            Sr = S_sum[1:] / np.maximum(counts[1:], 1.0)
            k_max = float(ks[int(np.argmax(Sr))])
            sel = ks <= 2.0 * k_max
            a_sf = float(np.trapezoid(Sr[sel], ks[sel]) / np.trapezoid(ks[sel] * Sr[sel], ks[sel]))
            out[target_steps[step]] = a_sf
    return out


def run(inputs, settings):
    M = float(inputs["M"])
    kappa = float(inputs["kappa"])
    gamma_dot = float(inputs["gamma_dot"])
    L = float(inputs["L"])
    phi_avg = float(inputs["phi_avg"])
    noise_amp = float(inputs["noise_amp"])
    t_targets = [float(t) for t in inputs["ref_a"]["x"]]

    flow_type = str(settings.get("flow_type", "shear"))
    N = int(settings.get("N", 128))
    dt = float(settings.get("dt", 0.05))
    n_ens = int(settings.get("n_ensemble", 3))
    seed = int(settings.get("seed", 0))

    acc = np.zeros(len(t_targets), dtype=np.float64)
    for j in range(n_ens):
        acc += _single_at(M, kappa, gamma_dot, L, phi_avg, noise_amp, flow_type, N, dt, seed + j, t_targets)
    acc /= float(n_ens)
    return {"a_sf": {"x": t_targets, "y": [float(v) for v in acc]}, "a_sf_final": float(acc[-1])}
