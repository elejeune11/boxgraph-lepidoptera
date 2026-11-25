import math
import numpy as np


def run(inputs, settings):
    chiA = float(inputs["chi_A"])
    chiB = float(inputs["chi_B"])
    gammaA = float(inputs["gamma_A"])
    kappa = float(inputs["kappa"])
    delta = float(inputs["delta"])
    L = float(inputs["L"])
    phi0A = float(inputs["phi0_A"])
    phi0B = float(inputs["phi0_B"])

    N = int(round(float(settings.get("N", 64))))
    dt = float(settings.get("dt", 2.0))
    t_max = float(settings.get("t_max", 20000.0))
    seed = int(settings.get("seed", 42))

    # Fourier-Galerkin pseudospectral solver (SI Eq. S4) with L-stable 2nd-order SBDF2 IMEX
    # and 2/3-rule dealiasing on the cubic nonlinearity.
    k = 2.0 * math.pi * np.fft.rfftfreq(N, d=L / N)
    k2 = k * k
    kAB = kappa - delta
    kBA = kappa + delta
    dealias = np.abs(np.fft.rfftfreq(N) * N) <= (N // 3)

    x = np.linspace(0.0, L, N, endpoint=False)
    rng = np.random.default_rng(seed)
    phiA = phi0A + 1.0e-2 * np.cos(2.0 * math.pi * x / L) + 1.0e-3 * rng.standard_normal(N)
    phiB = phi0B - 1.0e-2 * np.cos(2.0 * math.pi * x / L + 0.05) + 1.0e-3 * rng.standard_normal(N)
    phiA -= phiA.mean() - phi0A
    phiB -= phiB.mean() - phi0B

    Ahat = np.fft.rfft(phiA)
    Bhat = np.fft.rfft(phiB)

    M11 = -k2 * (chiA + gammaA * k2)
    M12 = -k2 * kAB
    M21 = -k2 * kBA
    M22 = -k2 * chiB

    # First step: L-stable backward-Euler IMEX (I - dt*M)
    e11 = 1.0 - dt * M11
    e12 = -dt * M12
    e21 = -dt * M21
    e22 = 1.0 - dt * M22
    edet = e11 * e22 - e12 * e21

    # Subsequent steps: L-stable SBDF2 (3*I - 2*dt*M) u^{n+1} = 4*u^n - u^{n-1} + 2*dt*(2*N^n - N^{n-1})
    s11 = 3.0 - 2.0 * dt * M11
    s12 = -2.0 * dt * M12
    s21 = -2.0 * dt * M21
    s22 = 3.0 - 2.0 * dt * M22
    sdet = s11 * s22 - s12 * s21

    steps = int(round(t_max / dt))
    sample_start = max(0, steps - int(round(2000.0 / dt)))
    sample_every = max(1, int(round(20.0 / dt)))

    A = np.fft.irfft(Ahat, n=N)
    nlA_prev = -(1.0 / 3.0) * k2 * np.fft.rfft(A * A * A) * dealias
    Ahat_prev = Ahat.copy()
    Bhat_prev = Bhat.copy()

    # Step 0 (IMEX1)
    rhsA = Ahat + dt * nlA_prev
    rhsB = Bhat
    Ahat = (e22 * rhsA - e12 * rhsB) / edet
    Bhat = (-e21 * rhsA + e11 * rhsB) / edet

    phases = []
    times = []
    for n in range(1, steps):
        A = np.fft.irfft(Ahat, n=N)
        nlA = -(1.0 / 3.0) * k2 * np.fft.rfft(A * A * A) * dealias

        rhsA = 4.0 * Ahat - Ahat_prev + 2.0 * dt * (2.0 * nlA - nlA_prev)
        rhsB = 4.0 * Bhat - Bhat_prev

        Ahat_prev = Ahat
        Bhat_prev = Bhat
        nlA_prev = nlA

        Ahat = (s22 * rhsA - s12 * rhsB) / sdet
        Bhat = (-s21 * rhsA + s11 * rhsB) / sdet

        if n >= sample_start and (n - sample_start) % sample_every == 0:
            phases.append(math.atan2(Ahat[1].imag, Ahat[1].real))
            times.append(n * dt)

    rhoA = float(abs(Ahat[1]) / N)
    rhoB = float(abs(Bhat[1]) / N)
    rhoA3 = float(abs(Ahat[3]) / N) if len(Ahat) > 3 else 0.0
    rhoB3 = float(abs(Bhat[3]) / N) if len(Bhat) > 3 else 0.0

    if rhoA < 1.0e-4 and rhoB < 1.0e-4:
        return {
            "rho_A": rhoA,
            "rho_B": rhoB,
            "rho_A_3": rhoA3,
            "rho_B_3": rhoB3,
            "v": 0.0,
            "psi": 0.0,
            "amp_ratio": None,
            "state_code": 0.0,
        }

    ph = np.unwrap(phases)
    q1 = 2.0 * math.pi / L
    v = float(abs((ph[-1] - ph[0]) / (times[-1] - times[0])) / q1)
    dth = math.atan2(Ahat[1].imag, Ahat[1].real) - math.atan2(Bhat[1].imag, Bhat[1].real)
    dth = (dth + math.pi) % (2.0 * math.pi) - math.pi
    psi = float(abs(abs(dth) - math.pi))
    state = 2.0 if v > 1.0e-4 else 1.0
    return {
        "rho_A": rhoA,
        "rho_B": rhoB,
        "rho_A_3": rhoA3,
        "rho_B_3": rhoB3,
        "v": v,
        "psi": psi,
        "amp_ratio": float(rhoA / rhoB) if rhoB > 1e-12 else None,
        "state_code": state,
    }
