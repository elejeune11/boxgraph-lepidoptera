import math
import numpy as np
from scipy.ndimage import correlate1d, label

_W8 = np.array(
    [1.0 / 280.0, -4.0 / 105.0, 1.0 / 5.0, -4.0 / 5.0, 0.0, 4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0],
    dtype=np.float64,
)


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    D = float(inputs["D"])
    phi_0 = float(inputs["phi_0"])

    L = int(settings.get("L", 32))
    dx = float(settings.get("dx", 1.0))
    dt = float(settings.get("dt", 0.001))
    t_end = float(settings.get("t_end", 2.0))
    seed = int(settings.get("seed", 0))

    rng = np.random.default_rng(seed)
    phi = np.full((L, L), phi_0, dtype=np.float64)
    noise_amp = math.sqrt(2.0 * D * dt / (dx * dx))
    nsteps = int(round(t_end / dt))

    for _ in range(nsteps):
        px = correlate1d(phi, _W8, axis=0, mode="wrap") / dx
        py = correlate1d(phi, _W8, axis=1, mode="wrap") / dx
        lap = correlate1d(px, _W8, axis=0, mode="wrap") / dx + correlate1d(py, _W8, axis=1, mode="wrap") / dx
        mu = A * (phi**3 - phi) - K * lap + lam * (px * px + py * py)
        Jx = -correlate1d(mu, _W8, axis=0, mode="wrap") / dx + zeta * lap * px
        Jy = -correlate1d(mu, _W8, axis=1, mode="wrap") / dx + zeta * lap * py
        divJ = correlate1d(Jx, _W8, axis=0, mode="wrap") / dx + correlate1d(Jy, _W8, axis=1, mode="wrap") / dx
        gx = rng.standard_normal((L, L))
        gy = rng.standard_normal((L, L))
        divG = correlate1d(gx, _W8, axis=0, mode="wrap") / dx + correlate1d(gy, _W8, axis=1, mode="wrap") / dx
        phi = phi - dt * divJ - noise_amp * divG
        if not np.isfinite(phi).all() or float(np.max(np.abs(phi))) > 1e3:
            raise ValueError("Verbatim Appendix A Eqs. (A1)-(A2) blew up in finite time")

    labeled, n_comp = label(phi < 0.0)
    radii = []
    for k in range(1, n_comp + 1):
        area = float(np.sum(labeled == k)) * (dx * dx)
        if area >= 4.0 * dx * dx:
            radii.append(math.sqrt(area / math.pi))
    R_bar = float(np.mean(radii)) if radii else 0.0
    n_bar = float(len(radii)) / ((L * dx) ** 2)
    return {
        "phi_min": float(np.min(phi)),
        "phi_max": float(np.max(phi)),
        "R_bar": R_bar,
        "n_bar": n_bar,
    }
