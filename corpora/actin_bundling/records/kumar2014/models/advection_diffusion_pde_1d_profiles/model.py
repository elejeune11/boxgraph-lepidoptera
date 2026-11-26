import math
import numpy as np


def _periodic_interp(x_eval, L, grid_vals):
    """Spectrally accurate Fourier interpolation of a periodic profile onto arbitrary x_eval."""
    N = len(grid_vals)
    coeffs = np.fft.rfft(grid_vals) / float(N)
    freqs = np.arange(len(coeffs), dtype=float)
    x_arr = np.asarray(x_eval, dtype=float)
    phase = (2.0 * math.pi / L) * np.outer(x_arr, freqs)
    weights = np.full(len(coeffs), 2.0, dtype=float)
    weights[0] = 1.0
    if N % 2 == 0:
        weights[-1] = 1.0
    out = np.real(np.exp(1j * phase) @ (coeffs * weights))
    return [float(v) for v in out]


def simulate_profiles(inputs, settings, report_unscaled: bool = False):
    P = float(inputs["P"])
    beta = float(inputs["beta"])
    alpha = float(inputs["alpha"])
    xi = float(inputs["xi"])
    f_0 = float(inputs["f_0"])
    L = float(inputs["L_over_l"])
    t_max = float(inputs["t_max"])
    eps = float(inputs["eps_init"])
    phi_0 = float(inputs["phi_0"])

    N = int(settings.get("n_grid", 64))
    dt = float(settings.get("dt", 0.03))

    dx = L / N
    x = np.arange(N, dtype=float) * dx
    k = 2.0 * math.pi * np.fft.rfftfreq(N, d=dx)
    k2 = k * k
    dealias = np.abs(k) <= (2.0 / 3.0) * (math.pi / dx)
    scale_f = ((1.0 + xi) ** 2) / xi
    ik = 1j * k
    vel_op = P * ik / (1.0 + k2)
    ik_d = -ik * dealias

    # Symmetric mode-1 initial perturbation centered at x = L/2 fixes the standing-wave spatial phase
    A = 1.0 + eps * np.cos(2.0 * math.pi * (x - 0.5 * L) / L)
    I = 1.0 + eps * np.cos(2.0 * math.pi * (x - 0.5 * L) / L)
    Ah = np.fft.rfft(A)
    Ih = np.fft.rfft(I)

    eA = np.exp(-k2 * dt)
    eI = np.exp(-alpha * k2 * dt)
    phi1A = np.ones_like(k2)
    phi1I = np.ones_like(k2)
    nz = k2 > 0.0
    phi1A[nz] = (1.0 - eA[nz]) / (k2[nz] * dt)
    phi1I[nz] = (1.0 - eI[nz]) / (alpha * k2[nz] * dt)

    n_steps = int(round(t_max / dt))
    t_record_start = max(0.0, t_max - 110.0)

    t_rec, A_rec, I_rec, v_rec, f_rec = [], [], [], [], []
    for step in range(n_steps + 1):
        Ar = np.fft.irfft(Ah, n=N)
        Ir = np.fft.irfft(Ih, n=N)
        fr_unsc = f_0 + (1.0 + beta) * Ar / (Ar + xi) + (1.0 - beta) * Ir / (Ir + xi)
        fr = scale_f * fr_unsc
        vr = np.fft.irfft(vel_op * np.fft.rfft(fr), n=N)
        nlA1 = ik_d * np.fft.rfft(vr * Ar)
        nlI1 = ik_d * np.fft.rfft(vr * Ir)

        t_cur = step * dt
        if t_cur >= t_record_start:
            t_rec.append(t_cur)
            A_rec.append(Ar)
            I_rec.append(Ir)
            v_rec.append(vr)
            f_rec.append(fr_unsc if report_unscaled else fr)
        if step == n_steps:
            break

        Ah_p = eA * Ah + dt * phi1A * nlA1
        Ih_p = eI * Ih + dt * phi1I * nlI1
        Ar2 = np.fft.irfft(Ah_p, n=N)
        Ir2 = np.fft.irfft(Ih_p, n=N)
        fr2 = scale_f * (f_0 + (1.0 + beta) * Ar2 / (Ar2 + xi) + (1.0 - beta) * Ir2 / (Ir2 + xi))
        vr2 = np.fft.irfft(vel_op * np.fft.rfft(fr2), n=N)
        nlA2 = ik_d * np.fft.rfft(vr2 * Ar2)
        nlI2 = ik_d * np.fft.rfft(vr2 * Ir2)

        Ah = eA * Ah + 0.5 * dt * (eA * nlA1 + nlA2)
        Ih = eI * Ih + 0.5 * dt * (eI * nlI1 + nlI2)

    t_rec = np.array(t_rec)
    A_rec = np.array(A_rec)
    I_rec = np.array(I_rec)
    v_rec = np.array(v_rec)
    f_rec = np.array(f_rec)

    I_mid = I_rec[:, N // 2]
    pks = np.where((I_mid[1:-1] > I_mid[:-2]) & (I_mid[1:-1] >= I_mid[2:]) & (I_mid[1:-1] > 2.0))[0] + 1
    if len(pks) < 2:
        # Fallback if a perturbation parameter pushes the system into a non-oscillatory state
        snap_idx = len(t_rec) - 1
        times_idx = {"a": snap_idx, "c": snap_idx, "e": snap_idx}
        res = {}
        for top_p, bot_p in [("a", "b"), ("c", "d"), ("e", "f")]:
            ii = times_idx[top_p]
            for out_key, ref_key, snap in [
                (f"A_{top_p}", f"A_{top_p}_ref", A_rec[ii]),
                (f"I_{top_p}", f"I_{top_p}_ref", I_rec[ii]),
                (f"v_{bot_p}", f"v_{bot_p}_ref", v_rec[ii]),
                (f"f_{bot_p}", f"f_{bot_p}_ref", f_rec[ii]),
            ]:
                xs = [float(xv) for xv in inputs[ref_key]["x"]]
                res[out_key] = {"x": xs, "y": _periodic_interp(xs, L, snap)}
        return res

    T_full = float(t_rec[pks[-1]] - t_rec[pks[-2]])
    t_pk = None
    for idx in reversed(pks):
        tp = float(t_rec[idx])
        if (tp + phi_0 * T_full >= t_rec[0]) and (tp + (phi_0 + 0.5) * T_full <= t_rec[-1]):
            t_pk = tp
            break
    if t_pk is None:
        t_pk = float(t_rec[pks[-2]])

    def interp_time(arr, t_target):
        j = int(np.searchsorted(t_rec, t_target))
        if j <= 0:
            return arr[0]
        if j >= len(t_rec):
            return arr[-1]
        w = (t_target - t_rec[j - 1]) / (t_rec[j] - t_rec[j - 1])
        return (1.0 - w) * arr[j - 1] + w * arr[j]

    times = {
        "a": t_pk + phi_0 * T_full,
        "c": t_pk + (phi_0 + 0.25) * T_full,
        "e": t_pk + (phi_0 + 0.50) * T_full,
    }
    res = {}
    for top_p, bot_p in [("a", "b"), ("c", "d"), ("e", "f")]:
        tt = times[top_p]
        A_snap = interp_time(A_rec, tt)
        I_snap = interp_time(I_rec, tt)
        v_snap = interp_time(v_rec, tt)
        f_snap = interp_time(f_rec, tt)
        for out_key, ref_key, snap in [
            (f"A_{top_p}", f"A_{top_p}_ref", A_snap),
            (f"I_{top_p}", f"I_{top_p}_ref", I_snap),
            (f"v_{bot_p}", f"v_{bot_p}_ref", v_snap),
            (f"f_{bot_p}", f"f_{bot_p}_ref", f_snap),
        ]:
            xs = [float(xv) for xv in inputs[ref_key]["x"]]
            res[out_key] = {"x": xs, "y": _periodic_interp(xs, L, snap)}
    return res


def run(inputs, settings):
    return simulate_profiles(inputs, settings, report_unscaled=False)
