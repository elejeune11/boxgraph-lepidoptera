import math
import numpy as np
from scipy.optimize import fsolve


def _flat_binodals(A, K, zeta, lam):
    alpha = (zeta - 2.0 * lam) / K
    nodes, weights = np.polynomial.legendre.leggauss(48)

    def pg(phi):
        h = 0.5 * phi
        xs = h * (nodes + 1.0)
        ex = np.exp(alpha * xs)
        return float(h * np.dot(weights, ex)), float(h * np.dot(weights, A * (xs**3 - xs) * ex))

    def res(x):
        p1, p2 = float(x[0]), float(x[1])
        m1, m2 = A * (p1**3 - p1), A * (p2**3 - p2)
        s1, g1 = pg(p1)
        s2, g2 = pg(p2)
        return [m1 - m2, (s1 * m1 - g1) - (s2 * m2 - g2)]

    p1, p2 = fsolve(res, [-1.0, 1.0], xtol=1e-12)
    return min(float(p1), float(p2)), max(float(p1), float(p2))


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    R1_0 = float(inputs["R1_0"])
    R2_0 = float(inputs["R2_0"])

    L = int(settings.get("L", 64))
    dx = float(settings.get("dx", 1.0))
    dt = float(settings.get("dt", 0.05))
    t_relax = float(settings.get("t_relax", 150.0))
    t_end = float(settings.get("t_end", 500.0))

    p1, p2 = _flat_binodals(A, K, zeta, lam)
    L_phys = L * dx
    x = np.arange(L, dtype=np.float64) * dx
    X, Y = np.meshgrid(x, x, indexing="ij")

    def dist(cx, cy):
        dx_ = (X - cx + 0.5 * L_phys) % L_phys - 0.5 * L_phys
        dy_ = (Y - cy + 0.5 * L_phys) % L_phys - 0.5 * L_phys
        return np.hypot(dx_, dy_)

    # Scale droplet radii proportionally if input radii exceed half-box clearance
    scale = min(1.0, (0.20 * L_phys) / max(R1_0, 1e-6))
    r1_init = R1_0 * scale
    r2_init = R2_0 * scale
    r1 = dist(0.25 * L_phys, 0.5 * L_phys)
    r2 = dist(0.75 * L_phys, 0.5 * L_phys)
    xi = math.sqrt(K / (2.0 * A))
    phi = (
        p1
        + 0.5 * (p2 - p1) * (1.0 - np.tanh((r1 - r1_init) / (1.5 * xi)))
        + 0.5 * (p2 - p1) * (1.0 - np.tanh((r2 - r2_init) / (1.5 * xi)))
    )

    # Exact Fourier symbols of the 8th-order central finite-difference operator Eq. (A2)
    kx_raw = 2.0 * np.pi * np.fft.fftfreq(L, d=dx)
    ky_raw = 2.0 * np.pi * np.fft.rfftfreq(L, d=dx)

    def keff(k):
        return (2.0 / dx) * (
            (4.0 / 5.0) * np.sin(k * dx)
            - (1.0 / 5.0) * np.sin(2.0 * k * dx)
            + (4.0 / 105.0) * np.sin(3.0 * k * dx)
            - (1.0 / 280.0) * np.sin(4.0 * k * dx)
        )

    ikx = 1j * keff(kx_raw)[:, None]
    iky = 1j * keff(ky_raw)[None, :]
    lap_hat = (ikx**2 + iky**2).real
    k2 = -lap_hat
    S_stab = 2.0 * A
    denom = 1.0 + dt * (K * k2**2 + S_stab * k2)

    def radii(p):
        mL = float(np.sum(np.clip((p[: L // 2, :] - p1) / (p2 - p1), 0.0, 2.0))) * (dx * dx)
        mR = float(np.sum(np.clip((p[L // 2 :, :] - p1) / (p2 - p1), 0.0, 2.0))) * (dx * dx)
        return math.sqrt(max(mL, 0.0) / math.pi), math.sqrt(max(mR, 0.0) / math.pi)

    phi_hat = np.fft.rfft2(phi)
    nsteps = int(round(t_end / dt))
    n_relax = int(round(t_relax / dt))
    R1_ref, R2_ref = radii(phi)

    for step in range(nsteps):
        px = np.fft.irfft2(ikx * phi_hat, s=(L, L))
        py = np.fft.irfft2(iky * phi_hat, s=(L, L))
        lap = np.fft.irfft2(lap_hat * phi_hat, s=(L, L))
        mu_nl_hat = np.fft.rfft2(A * (phi**3 - phi) + lam * (px * px + py * py))
        div_Jzeta_hat = ikx * np.fft.rfft2(zeta * lap * px) + iky * np.fft.rfft2(zeta * lap * py)
        phi_hat = (phi_hat + dt * (lap_hat * mu_nl_hat - div_Jzeta_hat + S_stab * k2 * phi_hat)) / denom
        phi = np.fft.irfft2(phi_hat, s=(L, L))
        if step == n_relax:
            R1_ref, R2_ref = radii(phi)

    R1_f, R2_f = radii(phi)
    px = np.fft.irfft2(ikx * phi_hat, s=(L, L))
    py = np.fft.irfft2(iky * phi_hat, s=(L, L))
    lap = np.fft.irfft2(lap_hat * phi_hat, s=(L, L))
    mu_hat = np.fft.rfft2(A * (phi**3 - phi) + lam * (px * px + py * py)) - K * lap_hat * phi_hat
    Jx = -np.fft.irfft2(ikx * mu_hat, s=(L, L)) + zeta * lap * px
    Jy = -np.fft.irfft2(iky * mu_hat, s=(L, L)) + zeta * lap * py

    Jx_mid = float(Jx[L // 2, L // 2])
    vapor_mask = phi < p1 + 0.2 * (p2 - p1)
    J_mag_vapor = float(np.median(np.hypot(Jx, Jy)[vapor_mask])) if np.any(vapor_mask) else float(abs(Jx_mid))
    dR1 = R1_f - R1_ref
    dR2 = R2_f - R2_ref
    # Forward Ostwald: larger droplet (left) grows relative to smaller droplet (right) -> (R1 - R2) increases (+1)
    ostwald_sign = 1.0 if (R1_f - R2_f) > (R1_ref - R2_ref) else -1.0
    return {
        "dR1": dR1,
        "dR2": dR2,
        "Jx_mid": Jx_mid,
        "J_mag_vapor": J_mag_vapor,
        "ostwald_sign": ostwald_sign,
    }
