import numpy as np


def run(inputs, settings):
    N = int(round(settings.get("N", 64)))
    n0 = int(round(inputs["n0"]))
    dphi = float(inputs["dphi"])

    x = np.arange(N, dtype=float) / N
    xx, _ = np.meshgrid(x, x, indexing="ij")
    I0 = np.cos(2.0 * np.pi * n0 * xx)
    It = np.cos(2.0 * np.pi * n0 * xx - dphi)
    dI = It - I0

    # Corrected DDM (matching Fig. S3 y-axis and Refs. 32-33): |FFT[I(x, t) - I(x, 0)]|^2
    D_2d = np.abs(np.fft.fft2(dI)) ** 2

    kx = np.fft.fftfreq(N) * N
    ky = np.fft.fftfreq(N) * N
    kxx, kyy = np.meshgrid(kx, ky, indexing="ij")
    n_rad = np.rint(np.hypot(kxx, kyy)).astype(int)

    mask_k0 = n_rad == n0
    mask_2k0 = n_rad == (2 * n0)

    D_k0_real = float(np.mean(D_2d[mask_k0]))
    D_2k0_real = float(np.mean(D_2d[mask_2k0]))
    D_k0_ratio = 1.0 if D_k0_real > 0.0 else 0.0

    return {
        "D_k0_real": D_k0_real,
        "D_k0_ratio": D_k0_ratio,
        "D_2k0_real": D_2k0_real,
    }
