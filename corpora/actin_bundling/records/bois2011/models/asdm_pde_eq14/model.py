import math
import numpy as np


def run(inputs, settings):
    Pe = float(inputs["Pe"])
    d = float(inputs["Da_over_Ds"])
    rho_a = float(inputs["rho_a"])
    rho_s = float(inputs["rho_s"])
    L = float(inputs["L_over_l"])
    T_end = float(inputs["T_end"])
    eps = float(inputs["eps_init"])
    N = int(round(settings.get("n_grid", 64)))
    dt = float(settings.get("dt", 0.005))

    X = np.linspace(0.0, L, N, endpoint=False)
    k = 2.0 * math.pi * np.fft.rfftfreq(N, d=L / N)
    # 3-mode periodic perturbation consistent with the 3-peak steady state in Fig. 5(d)
    a = 1.0 + eps * np.cos(3.0 * 2.0 * math.pi * X / L)
    s = 1.0 - eps * np.cos(3.0 * 2.0 * math.pi * X / L)
    a_hat = np.fft.rfft(a)
    s_hat = np.fft.rfft(s)

    nsteps = max(1, int(round(T_end / dt)))
    dt_actual = T_end / nsteps
    lin_a = np.exp(-d * (k**2) * dt_actual)
    lin_s = np.exp(-(k**2) * dt_actual)
    vel_op = 1j * k * Pe / (1.0 + k**2)
    cutoff = (len(k) * 2) // 3

    def rhs_nl(a_arr, s_arr):
        f_hat = np.fft.rfft(a_arr / (1.0 + a_arr))
        V = np.fft.irfft(vel_op * f_hat, n=N)
        a2s = a_arr * a_arr * s_arr
        Na = -1j * k * np.fft.rfft(V * a_arr) + np.fft.rfft(rho_a * (a2s - a_arr))
        Ns = -1j * k * np.fft.rfft(V * s_arr) + np.fft.rfft(rho_s * (1.0 - a2s))
        Na[cutoff:] = 0.0
        Ns[cutoff:] = 0.0
        return Na, Ns, V

    for _ in range(nsteps):
        a_cur = np.fft.irfft(a_hat, n=N)
        s_cur = np.fft.irfft(s_hat, n=N)
        Na1, Ns1, _ = rhs_nl(a_cur, s_cur)
        a_p = np.fft.irfft(lin_a * (a_hat + dt_actual * Na1), n=N)
        s_p = np.fft.irfft(lin_s * (s_hat + dt_actual * Ns1), n=N)
        Na2, Ns2, _ = rhs_nl(a_p, s_p)
        a_hat = lin_a * a_hat + 0.5 * dt_actual * (lin_a * Na1 + Na2)
        s_hat = lin_s * s_hat + 0.5 * dt_actual * (lin_s * Ns1 + Ns2)

    a_f = np.fft.irfft(a_hat, n=N)
    s_f = np.fft.irfft(s_hat, n=N)
    _, _, V_f = rhs_nl(a_f, s_f)
    v_over_U_max = float(np.max(np.abs(V_f)) / Pe) if Pe != 0.0 else 0.0
    peaks = int(np.sum((a_f > np.roll(a_f, 1)) & (a_f > np.roll(a_f, -1)) & (a_f > 2.0)))
    return {
        "a_max": float(np.max(a_f)),
        "a_min": float(np.min(a_f)),
        "s_max": float(np.max(s_f)),
        "s_min": float(np.min(s_f)),
        "v_over_U_max": v_over_U_max,
        "n_peaks": float(peaks),
    }
