import math
import numpy as np


def run(inputs, settings):
    beta = float(inputs["beta"])
    chi = float(inputs["chi"])
    gamma = float(inputs["gamma"])
    E1 = float(inputs["E1"])
    E3 = float(inputs["E3"])
    phi0 = float(inputs["phi0"])
    noise_var = float(inputs["noise_var"])
    k_cut_ratio = float(inputs["k_cut_ratio"])
    t_fit_start = float(inputs["t_fit_start"])
    t_end = float(inputs["t_end"])

    N = int(settings.get("n_grid", 64))
    L = float(settings.get("L_box", 62.5))
    dt = float(settings.get("dt", 0.4))
    seed = int(settings.get("seed", 0))
    n_samples = int(settings.get("n_samples", 6))

    rng = np.random.default_rng(seed)
    dx = L / N
    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    ky = 2.0 * np.pi * np.fft.rfftfreq(N, d=dx)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K2 = KX**2 + KY**2
    K4 = K2**2
    K6 = K2**3
    k_max_nyq = np.pi / dx
    dealias = (np.abs(KX) < (2.0 / 3.0) * k_max_nyq) & (np.abs(KY) < (2.0 / 3.0) * k_max_nyq)

    iKX = 1j * KX
    iKY = 1j * KY
    mKX2 = -KX**2
    mKY2 = -KY**2
    mKXKY = -KX * KY

    S = 1.0 / math.sqrt(beta)
    lam = 1.0 + 0.5 * abs(chi)
    sigma = math.sqrt(noise_var)

    phi = phi0 + rng.normal(0.0, sigma, (N, N))
    Qxx = rng.normal(0.0, sigma, (N, N))
    Qxy = rng.normal(0.0, sigma, (N, N))

    phi_k = np.fft.rfft2(phi)
    Qxx_k = np.fft.rfft2(Qxx)
    Qxy_k = np.fft.rfft2(Qxy)

    # IMEX factors from SI Eqs. (S31)-(S32) with phi^3 from Eqs. (9)/(S19)
    denom_phi = 1.0 / (1.0 + (gamma * lam * K4 + K6) * dt)
    num_phi_fac = 1.0 + K6 * dt
    denom_Q = 1.0 / (1.0 + E1 * K2 * dt)

    # Precompute radial binning tables for Eq. (S1) structure-factor domain size
    w_rfft = np.full(K2.shape, 2.0)
    w_rfft[:, 0] = 1.0
    if N % 2 == 0:
        w_rfft[:, -1] = 1.0
    k_mag = np.sqrt(K2)
    dk = 2.0 * np.pi / L
    nbins = N // 2
    bin_idx = np.rint(k_mag / dk).astype(int)
    valid = (bin_idx >= 1) & (bin_idx < nbins)
    idx_valid = bin_idx[valid]
    w_valid = w_rfft[valid]
    cnt = np.bincount(idx_valid, weights=w_valid, minlength=nbins)
    nz = cnt > 0
    k_r = np.arange(nbins) * dk

    def compute_R(pk_field):
        Pk = np.abs(pk_field) ** 2
        S_r = np.bincount(idx_valid, weights=(Pk * w_rfft)[valid], minlength=nbins)
        S_avg = np.zeros(nbins)
        S_avg[nz] = S_r[nz] / cnt[nz]
        k_peak_idx = int(np.argmax(S_avg[1:])) + 1
        k_cut = k_cut_ratio * k_r[k_peak_idx]
        sel = (k_r > 0) & (k_r <= k_cut)
        den = float(np.sum(k_r[sel] * S_avg[sel]))
        return float(np.sum(S_avg[sel]) / den) if den > 0.0 else float("nan")

    nsteps = int(round(t_end / dt))
    sample_times = np.geomspace(t_fit_start, t_end, n_samples)
    sample_steps = {int(round(ts / dt)): float(ts) for ts in sample_times}
    sample_steps[nsteps] = t_end

    ts_rec = []
    Rs_rec = []

    for step in range(1, nsteps + 1):
        dphi_x = np.fft.irfft2(iKX * phi_k, s=(N, N))
        dphi_y = np.fft.irfft2(iKY * phi_k, s=(N, N))
        trQ2 = 2.0 * (Qxx * Qxx + Qxy * Qxy)

        fx = Qxx * dphi_x + Qxy * dphi_y
        fy = Qxy * dphi_x - Qxx * dphi_y
        div_Qgrad = np.fft.irfft2((iKX * np.fft.rfft2(fx) + iKY * np.fft.rfft2(fy)) * dealias, s=(N, N))

        nl_phi = phi**3 - phi - 0.5 * trQ2 - (chi / S) * div_Qgrad
        phi_k_new = (phi_k * num_phi_fac - gamma * K2 * dt * (np.fft.rfft2(nl_phi) * dealias)) * denom_phi

        bulk_fac = 1.0 + phi - beta * trQ2
        grad2 = dphi_x * dphi_x + dphi_y * dphi_y
        anch_xx = -(chi / (2.0 * S)) * (dphi_x * dphi_x - 0.5 * grad2)
        anch_xy = -(chi / (2.0 * S)) * (dphi_x * dphi_y)

        if E3 != 0.0:
            dQxx_x = np.fft.irfft2(iKX * Qxx_k, s=(N, N))
            dQxx_y = np.fft.irfft2(iKY * Qxx_k, s=(N, N))
            dQxy_x = np.fft.irfft2(iKX * Qxy_k, s=(N, N))
            dQxy_y = np.fft.irfft2(iKY * Qxy_k, s=(N, N))
            d2Qxx_xx = np.fft.irfft2(mKX2 * Qxx_k, s=(N, N))
            d2Qxx_yy = np.fft.irfft2(mKY2 * Qxx_k, s=(N, N))
            d2Qxx_xy = np.fft.irfft2(mKXKY * Qxx_k, s=(N, N))
            d2Qxy_xx = np.fft.irfft2(mKX2 * Qxy_k, s=(N, N))
            d2Qxy_yy = np.fft.irfft2(mKY2 * Qxy_k, s=(N, N))
            d2Qxy_xy = np.fft.irfft2(mKXKY * Qxy_k, s=(N, N))

            # Printed Eq. (10) / Eq. (S20) cubic elastic term (+E_3/2 * (...))
            e3_xx = -E3 * (
                Qxx * (d2Qxx_xx - d2Qxx_yy)
                + 2.0 * Qxy * d2Qxx_xy
                + 0.5 * dQxx_x**2
                + dQxx_x * dQxy_y
                - 0.5 * dQxx_y**2
                + dQxx_y * dQxy_x
                - 0.5 * dQxy_x**2
                + 0.5 * dQxy_y**2
            )
            e3_xy = -E3 * (
                Qxx * (d2Qxy_xx - d2Qxy_yy)
                + 2.0 * Qxy * d2Qxy_xy
                - dQxx_x * dQxx_y
                + dQxx_x * dQxy_x
                - dQxx_y * dQxy_y
                + dQxy_x * dQxy_y
            )
        else:
            e3_xx = 0.0
            e3_xy = 0.0

        Qxx_k_new = (Qxx_k + dt * (np.fft.rfft2(bulk_fac * Qxx + anch_xx + e3_xx) * dealias)) * denom_Q
        Qxy_k_new = (Qxy_k + dt * (np.fft.rfft2(bulk_fac * Qxy + anch_xy + e3_xy) * dealias)) * denom_Q

        phi_k, Qxx_k, Qxy_k = phi_k_new, Qxx_k_new, Qxy_k_new
        phi = np.fft.irfft2(phi_k, s=(N, N))
        Qxx = np.fft.irfft2(Qxx_k, s=(N, N))
        Qxy = np.fft.irfft2(Qxy_k, s=(N, N))

        if not np.isfinite(phi).all() or float(np.max(np.abs(phi))) > 1e3:
            return {
                "alpha_fit": float("nan"),
                "R_final": float("nan"),
                "R_growth_ratio": float("nan"),
                "Q_max_final": float("nan"),
            }

        if step in sample_steps:
            ts_rec.append(step * dt)
            Rs_rec.append(compute_R(phi_k))

    ts_arr = np.array(ts_rec, dtype=float)
    Rs_arr = np.array(Rs_rec, dtype=float)
    alpha_fit = float(np.polyfit(np.log(ts_arr), np.log(Rs_arr), 1)[0])
    R_final = float(Rs_arr[-1])
    R_growth_ratio = float(Rs_arr[-1] / Rs_arr[0]) if Rs_arr[0] > 0 else float("nan")
    Q_max_final = float(np.max(np.hypot(Qxx, Qxy)))
    return {
        "alpha_fit": alpha_fit,
        "R_final": R_final,
        "R_growth_ratio": R_growth_ratio,
        "Q_max_final": Q_max_final,
    }
