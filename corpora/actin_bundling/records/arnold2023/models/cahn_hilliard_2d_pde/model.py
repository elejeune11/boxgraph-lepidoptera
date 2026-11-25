import numpy as np
from scipy import fft as spfft


def run(inputs, settings):
    M = float(inputs["M"])
    kappa = float(inputs["kappa"])
    gamma_dot = float(inputs["gamma_dot"])
    L = float(inputs["L"])
    phi_avg = float(inputs["phi_avg"])
    noise_amp = float(inputs["noise_amp"])

    flow_type = str(settings.get("flow_type", "passive"))
    N = int(settings.get("N", 128))
    dt = float(settings.get("dt", 0.05))
    t_start_fit = float(settings.get("t_start_fit", 20.0))
    t_end = float(settings.get("t_end", 50.0))
    sample_dt = float(settings.get("sample_dt", 2.5))
    seed = int(settings.get("seed", 0))

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
    elif flow_type == "extensional":
        vx = gamma_dot * Y
        vy = gamma_dot * X
    else:
        raise ValueError(f"unknown flow_type: {flow_type!r}")

    phi_k = spfft.rfft2(phi)
    denom = 1.0 + dt * M * kappa * K4

    # Precompute radial shell binning for static structure factor Eq. (47)
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

    nsteps = int(round(t_end / dt))
    sample_every = max(1, int(round(sample_dt / dt)))

    ts = []
    a_sfs = []
    a_aps = []
    k_max_last = float("nan")
    phi_cur = phi

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

        t_cur = step * dt
        if t_cur >= t_start_fit - 1.0e-9 and (step % sample_every == 0 or step == nsteps):
            phi_cur = spfft.irfft2(phi_k, s=(N, N))
            S2d = np.abs(phi_k) ** 2
            S_sum = np.bincount(bins_m, weights=(S2d[mask] * w_m), minlength=max_b + 1)
            Sr = S_sum[1:] / np.maximum(counts[1:], 1.0)
            k_max_last = float(ks[int(np.argmax(Sr))])
            sel = ks <= 2.0 * k_max_last
            if int(np.sum(sel)) >= 2:
                a_sf = float(np.trapezoid(Sr[sel], ks[sel]) / np.trapezoid(ks[sel] * Sr[sel], ks[sel]))
            else:
                a_sf = 1.0 / k_max_last

            bw = phi_cur < 0.0
            perim_edges = int(np.sum(bw != np.roll(bw, 1, axis=1))) + int(np.sum(bw != np.roll(bw, 1, axis=0)))
            a_ap = float((bw.sum() * dx * dx) / (max(perim_edges, 1) * dx))

            ts.append(t_cur)
            a_sfs.append(a_sf)
            a_aps.append(a_ap)

    ts_arr = np.asarray(ts, dtype=np.float64)
    alpha_sf = float(np.polyfit(np.log(ts_arr), np.log(np.asarray(a_sfs, dtype=np.float64)), 1)[0])
    alpha_ap = float(np.polyfit(np.log(ts_arr), np.log(np.asarray(a_aps, dtype=np.float64)), 1)[0])
    mean_phi_final = float(np.mean(spfft.irfft2(phi_k, s=(N, N))))
    return {
        "alpha_sf": alpha_sf,
        "alpha_ap": alpha_ap,
        "mean_phi_final": mean_phi_final,
        "k_max_final": k_max_last,
    }
