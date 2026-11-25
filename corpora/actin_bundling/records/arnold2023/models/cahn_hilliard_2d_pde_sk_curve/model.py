"""2D pseudospectral IMEX Cahn-Hilliard solver returning radially averaged <S(k)> on the reference k-grid."""

import numpy as np
from scipy import fft as spfft


def _single_sk(M, kappa, gamma_dot, L, phi_avg, noise_amp, flow_type, N, dt, t_eval, seed, k_targets):
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

    phi_k = spfft.rfft2(phi)
    denom = 1.0 + dt * M * kappa * K4
    nsteps = int(round(t_eval / dt))

    for _ in range(nsteps):
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

    # Continuous Parseval-normalized static structure factor S(k) = (L^2 / N^4) |phi_k|^2
    S2d = (np.abs(phi_k) ** 2) * ((L * L) / (float(N) ** 4))
    w_rfft = np.ones_like(KX)
    w_rfft[:, 1:-1] = 2.0
    K_mag = np.sqrt(K2)

    # Bin onto dk = 0.1 shells matching Fig. 3(c) discrete k spacing
    sk_vals = []
    for kc in k_targets:
        m = (K_mag >= kc - 0.05) & (K_mag < kc + 0.05)
        wsum = float(np.sum(w_rfft[m]))
        sk_vals.append(float(np.sum(S2d[m] * w_rfft[m]) / wsum) if wsum > 0.0 else 0.0)
    return np.asarray(sk_vals, dtype=np.float64)


def run(inputs, settings):
    M = float(inputs["M"])
    kappa = float(inputs["kappa"])
    gamma_dot = float(inputs["gamma_dot"])
    L = float(inputs["L"])
    phi_avg = float(inputs["phi_avg"])
    noise_amp = float(inputs["noise_amp"])
    k_targets = [float(k) for k in inputs["ref_Sk"]["x"]]

    flow_type = str(settings.get("flow_type", "shear"))
    N = int(settings.get("N", 128))
    dt = float(settings.get("dt", 0.05))
    t_eval = float(settings.get("t_eval", 50.0))
    n_ens = int(settings.get("n_ensemble", 3))
    seed = int(settings.get("seed", 0))

    acc = np.zeros(len(k_targets), dtype=np.float64)
    for j in range(n_ens):
        acc += _single_sk(M, kappa, gamma_dot, L, phi_avg, noise_amp, flow_type, N, dt, t_eval, seed + j, k_targets)
    acc /= float(n_ens)
    k_max = float(k_targets[int(np.argmax(acc))])
    return {"Sk": {"x": k_targets, "y": [float(v) for v in acc]}, "k_max": k_max}
