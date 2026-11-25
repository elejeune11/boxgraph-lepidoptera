import numpy as np


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])

    L = int(settings.get("L", 64))
    dx = float(settings.get("dx", 1.0))
    dt = float(settings.get("dt", 0.25))
    t_end = float(settings.get("t_end", 2000.0))

    # Periodic domain initialized with a dense slab on [L/4, 3L/4]
    x = np.arange(L, dtype=np.float64) * dx
    L_phys = L * dx
    phi = np.tanh((0.25 * L_phys - np.abs(x - 0.5 * L_phys)) / 2.0)

    # Exact Fourier symbol of the 8th-order central finite-difference operator Eq. (A2)
    k_raw = 2.0 * np.pi * np.fft.rfftfreq(L, d=dx)
    keff = (2.0 / dx) * (
        (4.0 / 5.0) * np.sin(k_raw * dx)
        - (1.0 / 5.0) * np.sin(2.0 * k_raw * dx)
        + (4.0 / 105.0) * np.sin(3.0 * k_raw * dx)
        - (1.0 / 280.0) * np.sin(4.0 * k_raw * dx)
    )
    ik = 1j * keff
    lap_hat = -(keff**2)
    k2 = keff**2

    # Semi-implicit (IMEX) splitting of linear -K d_x^4 phi + stabilization S d_x^2 phi;
    # uses the exact Eq. (A2) spatial operator and converges to Eq. (A1) as dt -> 0.
    S_stab = 2.0 * A
    denom = 1.0 + dt * (K * k2**2 + S_stab * k2)
    phi_hat = np.fft.rfft(phi)
    nsteps = int(round(t_end / dt))
    for _ in range(nsteps):
        px = np.fft.irfft(ik * phi_hat, n=L)
        lap = np.fft.irfft(lap_hat * phi_hat, n=L)
        mu_nl_hat = np.fft.rfft(A * (phi**3 - phi) + lam * (px * px))
        div_Jzeta_hat = ik * np.fft.rfft(zeta * lap * px)
        phi_hat = (phi_hat + dt * (lap_hat * mu_nl_hat - div_Jzeta_hat + S_stab * k2 * phi_hat)) / denom
        phi = np.fft.irfft(phi_hat, n=L)

    return {
        "phi_1": float(np.min(phi)),
        "phi_2": float(np.max(phi)),
    }
