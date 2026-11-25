"""Zero-padded 1D FFT dominant spatial frequency and ridge spacing (SI p. 11, Fig. S10D).

Resamples h(x) onto a uniform lateral grid (if nonuniform), subtracts mean, computes a
zero-padded real FFT of length n_fft, and selects the peak frequency f >= f_min (in 1/um).
"""

import numpy as np


def run(inputs, settings):
    curve = inputs["h"]
    x = np.asarray(curve["x"], dtype=float)
    y = np.asarray(curve["y"], dtype=float)
    if x.size < 4:
        return {"freq": float("nan"), "s": float("nan")}

    order = np.argsort(x)
    x, y = x[order], y[order]

    n_uniform = int(settings.get("n_uniform", 256))
    dx_arr = np.diff(x)
    if dx_arr.size == 0 or (np.max(dx_arr) - np.min(dx_arr)) > 1e-6 * np.mean(dx_arr):
        x_grid = np.linspace(x[0], x[-1], max(n_uniform, x.size))
        y_grid = np.interp(x_grid, x, y)
    else:
        x_grid, y_grid = x, y

    dx = float(x_grid[1] - x_grid[0])
    y_centered = y_grid - float(np.mean(y_grid))
    n_fft = max(int(settings.get("n_fft", 4096)), y_centered.size)
    spec = np.abs(np.fft.rfft(y_centered, n=n_fft))
    freqs = np.fft.rfftfreq(n_fft, d=dx)

    f_min = float(settings.get("f_min", 0.3))
    mask = freqs >= f_min
    if not np.any(mask):
        return {"freq": float("nan"), "s": float("nan")}

    sub_freqs = freqs[mask]
    sub_spec = spec[mask]
    peak_idx = int(np.argmax(sub_spec))
    freq_peak = float(sub_freqs[peak_idx])
    spacing = 1.0 / freq_peak if freq_peak > 0.0 else float("nan")
    return {"freq": freq_peak, "s": spacing}
