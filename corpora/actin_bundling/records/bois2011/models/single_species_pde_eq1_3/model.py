import math
import numpy as np


def run(inputs, settings):
    Pe = float(inputs["Pe"])
    L = float(inputs["L_over_l"])
    c0 = float(inputs["c0"])
    T_end = float(inputs["T_end"])
    eps = float(inputs["eps_init"])
    N = int(round(settings.get("n_grid", 64)))
    dt = float(settings.get("dt", 0.002))

    X = np.linspace(0.0, L, N, endpoint=False)
    k = 2.0 * math.pi * np.fft.rfftfreq(N, d=L / N)
    # Initial perturbation exciting both mode 2 (fastest-growing for L/l = 3*pi) and symmetry-breaking mode 1
    c = c0 * (
        1.0
        + eps * np.cos(2.0 * 2.0 * math.pi * X / L)
        + 0.1 * eps * np.cos(2.0 * math.pi * X / L + 0.5)
    )
    c_hat = np.fft.rfft(c)
    lin_full = np.exp(-(k**2) * dt)
    vel_op = 1j * k * Pe / (1.0 + k**2)
    cutoff = (len(k) * 2) // 3

    def adv_rhs(c_arr):
        f_hat = np.fft.rfft(c_arr / (1.0 + c_arr))
        V = np.fft.irfft(vel_op * f_hat, n=N)
        flux_hat = np.fft.rfft(V * c_arr)
        flux_hat[cutoff:] = 0.0
        return -1j * k * flux_hat

    nsteps = max(1, int(round(T_end / dt)))
    dt_actual = T_end / nsteps
    lin_full = np.exp(-(k**2) * dt_actual)

    for _ in range(nsteps):
        c_cur = np.fft.irfft(c_hat, n=N)
        N1 = adv_rhs(c_cur)
        c_pred = np.fft.irfft(lin_full * (c_hat + dt_actual * N1), n=N)
        N2 = adv_rhs(c_pred)
        c_hat = lin_full * c_hat + 0.5 * dt_actual * (lin_full * N1 + N2)

    c_final = np.fft.irfft(c_hat, n=N)
    # Count distinct concentration peaks significantly above c0
    thresh = 1.5 * c0
    peaks = int(
        np.sum((c_final > np.roll(c_final, 1)) & (c_final > np.roll(c_final, -1)) & (c_final > thresh))
    )
    return {
        "n_peaks_final": float(peaks),
        "c_max_final": float(np.max(c_final)),
        "c_mean_final": float(np.mean(c_final)),
    }
