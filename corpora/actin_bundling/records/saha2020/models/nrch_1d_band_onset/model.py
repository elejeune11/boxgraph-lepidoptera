import numpy as np


def _speed_at_alpha(alpha, c11, c12, c21, c22, chi, chip, kappa, p1b, p2b, N, L, dt, T, S):
    x = np.linspace(-0.5 * L, 0.5 * L, N, endpoint=False)
    a1 = abs(c12 - c11) / 2.0
    a2 = abs(c22 - c21) / 2.0
    # Single-domain bulk phase-separated profile with tiny parity-breaking seed
    p1 = p1b + a1 * np.tanh(12.0 * np.cos(2.0 * np.pi * x / L))
    p2 = p2b + a2 * np.tanh(12.0 * np.cos(2.0 * np.pi * x / L)) + 0.005 * np.sin(2.0 * np.pi * x / L)
    p1 += p1b - p1.mean()
    p2 += p2b - p2.mean()

    q = 2.0 * np.pi * np.fft.rfftfreq(N, d=L / N)
    q2 = q * q
    q4 = q2 * q2

    a11 = -S * q2 - kappa * q4
    a22 = -S * q2 - kappa * q4
    a12 = -(chi + alpha) * q2
    a21 = -(chi - alpha) * q2
    m11 = 1.0 - dt * a11
    m12 = -dt * a12
    m21 = -dt * a21
    m22 = 1.0 - dt * a22
    detM = m11 * m22 - m12 * m21
    inv11 = m22 / detM
    inv12 = -m12 / detM
    inv21 = -m21 / detM
    inv22 = m11 / detM

    p1h = np.fft.rfft(p1)
    p2h = np.fft.rfft(p2)
    nsteps = max(1, int(round(T / dt)))
    for _ in range(nsteps):
        r1 = 2.0 * (p1 - c11) * (p1 - c12) * (2.0 * p1 - c11 - c12) + 2.0 * chip * p1 * p2 * p2 - S * p1
        r2 = 2.0 * (p2 - c21) * (p2 - c22) * (2.0 * p2 - c21 - c22) + 2.0 * chip * p2 * p1 * p1 - S * p2
        rhs1 = p1h - dt * q2 * np.fft.rfft(r1)
        rhs2 = p2h - dt * q2 * np.fft.rfft(r2)
        p1h = inv11 * rhs1 + inv12 * rhs2
        p2h = inv21 * rhs1 + inv22 * rhs2
        p1 = np.fft.irfft(p1h, n=N)
        p2 = np.fft.irfft(p2h, n=N)

    dp1 = np.fft.irfft(1j * q * p1h, n=N)
    mu1_bulk = (
        2.0 * (p1 - c11) * (p1 - c12) * (2.0 * p1 - c11 - c12)
        + (chi + alpha) * p2
        + 2.0 * chip * p1 * p2 * p2
    )
    dp1_dt = np.fft.irfft(-q2 * np.fft.rfft(mu1_bulk) - kappa * q4 * p1h, n=N)
    denom = float(np.dot(dp1, dp1))
    return abs(float(np.dot(dp1_dt, dp1)) / denom) if denom > 1e-12 else 0.0


def run(inputs, settings):
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    chi = float(inputs["chi"])
    chip = float(inputs["chip"])
    kappa = float(inputs["kappa"])
    p1b = float(inputs["phi1_bar"])
    p2b = float(inputs["phi2_bar"])

    N = int(settings.get("N", 128))
    L = float(settings.get("L", 2.0))
    dt = float(settings.get("dt", 1e-3))
    T = float(settings.get("T", 12.0))
    S = float(settings.get("S", 0.8))
    v_thresh = float(settings.get("v_thresh", 0.05))
    n_bisect = int(settings.get("n_bisect", 8))

    lo = abs(chi)
    hi = abs(chi) + 0.06
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        v = _speed_at_alpha(mid, c11, c12, c21, c22, chi, chip, kappa, p1b, p2b, N, L, dt, T, S)
        if v > v_thresh:
            hi = mid
        else:
            lo = mid
    return {"alpha_onset": 0.5 * (lo + hi)}
