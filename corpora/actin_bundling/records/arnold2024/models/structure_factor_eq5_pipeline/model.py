import numpy as np


def run(inputs, settings):
    L = float(inputs["L"])
    n0 = int(round(inputs["n0"]))
    N = int(round(settings.get("N", 64)))
    min_wl_px = float(settings.get("min_wavelength_px", 5.0))

    # Deterministic 2D square pattern with spatial mode n0 along x and y (plus DC offset)
    x = np.arange(N, dtype=float) / N
    xx, yy = np.meshgrid(x, x, indexing="ij")
    img = 1.0 + np.cos(2.0 * np.pi * n0 * xx) + np.cos(2.0 * np.pi * n0 * yy)

    # As printed on p. 4: S(n) = |FFT(n)| (unsquared magnitude)
    S_2d = np.abs(np.fft.fft2(img))

    kx_idx = np.fft.fftfreq(N) * N
    ky_idx = np.fft.fftfreq(N) * N
    kxx, kyy = np.meshgrid(kx_idx, ky_idx, indexing="ij")
    n_rad = np.floor(np.hypot(kxx, kyy) + 0.5).astype(int)

    # Discard n = 0 and wavelengths 2*pi/k < min_wavelength_px pixels (i.e., N / n < min_wl_px)
    n_max = int(np.floor(N / min_wl_px))
    n_vals = np.arange(1, n_max + 1, dtype=int)
    S_radial = np.array([float(np.mean(S_2d[n_rad == n])) for n in n_vals], dtype=float)
    k_vals = 2.0 * np.pi * n_vals / L

    # Trapezoidal integration of Eq. (5): <k> = int k S(k) dk / int S(k) dk
    dk = np.diff(k_vals)
    y_num = k_vals * S_radial
    num = float(0.5 * np.sum((y_num[:-1] + y_num[1:]) * dk))
    den = float(0.5 * np.sum((S_radial[:-1] + S_radial[1:]) * dk))

    k_mean = num / den
    wavelength = 2.0 * np.pi / k_mean
    return {"k_mean": k_mean, "wavelength": wavelength}
