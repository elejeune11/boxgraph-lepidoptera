import math
import numpy as np
from scipy.ndimage import label


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    D = float(inputs["D"])
    phi_0 = float(inputs["phi_0"])

    L = int(settings.get("L", 48))
    dx = float(settings.get("dx", 1.0))
    dt = float(settings.get("dt", 0.05))
    t_end = float(settings.get("t_end", 200.0))
    k_cut = float(settings.get("k_cut", 1.0))
    seed = int(settings.get("seed", 0))

    rng = np.random.default_rng(seed)
    phi = np.full((L, L), phi_0, dtype=np.float64)
    kx = (2.0 * np.pi * np.fft.fftfreq(L, d=dx))[:, None]
    ky = (2.0 * np.pi * np.fft.rfftfreq(L, d=dx))[None, :]
    mask = np.hypot(kx, ky) <= k_cut
    ikx = 1j * kx * mask
    iky = 1j * ky * mask
    lap_hat = -(kx**2 + ky**2) * mask
    k2 = -lap_hat
    S_stab = 2.0 * A
    denom = 1.0 + dt * (K * k2**2 + S_stab * k2)
    noise_amp = math.sqrt(2.0 * D * dt / (dx * dx))

    phi_hat = np.fft.rfft2(phi)
    nsteps = int(round(t_end / dt))
    for _ in range(nsteps):
        px = np.fft.irfft2(ikx * phi_hat, s=(L, L))
        py = np.fft.irfft2(iky * phi_hat, s=(L, L))
        lap = np.fft.irfft2(lap_hat * phi_hat, s=(L, L))
        mu_nl_hat = np.fft.rfft2(A * (phi**3 - phi) + lam * (px * px + py * py)) * mask
        div_Jzeta_hat = ikx * np.fft.rfft2(zeta * lap * px) + iky * np.fft.rfft2(zeta * lap * py)
        gx = rng.standard_normal((L, L))
        gy = rng.standard_normal((L, L))
        div_G_hat = ikx * np.fft.rfft2(gx) + iky * np.fft.rfft2(gy)
        phi_hat = (
            phi_hat
            + dt * (lap_hat * mu_nl_hat - div_Jzeta_hat + S_stab * k2 * phi_hat)
            - noise_amp * div_G_hat
        ) / denom
        phi = np.fft.irfft2(phi_hat, s=(L, L))

    labeled, n_comp = label(phi < 0.0)
    radii = []
    for k in range(1, n_comp + 1):
        area = float(np.sum(labeled == k)) * (dx * dx)
        if area >= 4.0 * dx * dx:
            radii.append(math.sqrt(area / math.pi))
    R_bar = float(np.mean(radii)) if radii else 0.0
    n_bar = float(len(radii)) / ((L * dx) ** 2)
    return {
        "R_bar": R_bar,
        "n_bar": n_bar,
        "phi_min": float(np.min(phi)),
        "phi_max": float(np.max(phi)),
    }
