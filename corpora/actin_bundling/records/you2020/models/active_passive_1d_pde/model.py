import math
import numpy as np


def run(inputs, settings):
    v0 = float(inputs["v_0"])
    DR = float(inputs["D_R"])
    a0AA = float(inputs["a0_AA"])
    a0AP = float(inputs["a0_AP"])
    a1 = float(inputs["a1"])
    kap = float(inputs["kappa_ap"])
    rA0 = float(inputs["rho_A0"])
    rP0 = float(inputs["rho_P0"])
    L = float(inputs["L"])

    N = int(round(float(settings.get("N", 64))))
    dt = float(settings.get("dt", 0.05))
    t_max = float(settings.get("t_max", 800.0))
    seed = int(settings.get("seed", 1))

    k = 2.0 * math.pi * np.fft.rfftfreq(N, d=L / N)
    ik = 1j * k
    k4 = k**4
    dealias = np.abs(np.fft.rfftfreq(N) * N) <= (N // 3)

    x = np.linspace(0.0, L, N, endpoint=False)
    rng = np.random.default_rng(seed)
    rhoA = rA0 + 1.0e-2 * np.cos(2.0 * math.pi * x / L) + 1.0e-3 * rng.standard_normal(N)
    rhoP = rP0 - 1.0e-2 * np.cos(2.0 * math.pi * x / L + 0.05) + 1.0e-3 * rng.standard_normal(N)
    rhoA -= rhoA.mean() - rA0
    rhoP -= rhoP.mean() - rP0

    rA_hat = np.fft.rfft(rhoA)
    rP_hat = np.fft.rfft(rhoP)

    # L-stable SBDF2 on biharmonic -kappa*k^4 with 2/3-dealiased nonlinear flux
    e_lhs = 1.0 + dt * kap * k4
    s_lhs = 3.0 + 2.0 * dt * kap * k4

    def nl_flux(Ah, Ph):
        rA = np.fft.irfft(Ah, n=N)
        rP = np.fft.irfft(Ph, n=N)
        drA = np.fft.irfft(ik * Ah, n=N)
        drP = np.fft.irfft(ik * Ph, n=N)
        v = v0 * (1.0 - a0AA * rA - a0AP * rP)
        DAA = a1 * rA + (v / (2.0 * DR)) * (v - a0AA * v0 * rA)
        DAP = a1 * rA - (v / (2.0 * DR)) * a0AP * v0 * rA
        DPA = a1 * rP
        DPP = a1 * rP
        nA = ik * np.fft.rfft(DAA * drA + DAP * drP) * dealias
        nP = ik * np.fft.rfft(DPA * drA + DPP * drP) * dealias
        return nA, nP

    steps = int(round(t_max / dt))
    window_time = min(200.0, 0.3 * t_max)
    sample_start = max(0, steps - int(round(window_time / dt)))
    sample_every = max(1, int(round(2.0 / dt)))

    nA_prev, nP_prev = nl_flux(rA_hat, rP_hat)
    rA_prev = rA_hat.copy()
    rP_prev = rP_hat.copy()

    # Step 0 (IMEX1)
    rA_hat = (rA_hat + dt * nA_prev) / e_lhs
    rP_hat = (rP_hat + dt * nP_prev) / e_lhs

    phases = []
    times = []
    for n in range(1, steps):
        nA, nP = nl_flux(rA_hat, rP_hat)
        rhsA = 4.0 * rA_hat - rA_prev + 2.0 * dt * (2.0 * nA - nA_prev)
        rhsP = 4.0 * rP_hat - rP_prev + 2.0 * dt * (2.0 * nP - nP_prev)
        rA_prev = rA_hat
        rP_prev = rP_hat
        nA_prev = nA
        nP_prev = nP

        rA_hat = rhsA / s_lhs
        rP_hat = rhsP / s_lhs

        if n >= sample_start and (n - sample_start) % sample_every == 0:
            phases.append(math.atan2(rA_hat[1].imag, rA_hat[1].real))
            times.append(n * dt)

    rA_final = np.fft.irfft(rA_hat, n=N)
    amp_A = float(abs(rA_hat[1]) / N)
    if amp_A < 1.0e-3:
        return {
            "amp_A": amp_A,
            "rho_A_max": float(rA_final.max()),
            "rho_A_min": float(rA_final.min()),
            "omega_travel": 0.0,
            "period_travel": None,
            "state_code": 0.0,
        }

    ph = np.unwrap(phases)
    omega = float(abs((ph[-1] - ph[0]) / (times[-1] - times[0])))
    if omega > 1.0e-3:
        state = 2.0
        period = float(2.0 * math.pi / omega)
    else:
        state = 1.0
        period = None

    return {
        "amp_A": amp_A,
        "rho_A_max": float(rA_final.max()),
        "rho_A_min": float(rA_final.min()),
        "omega_travel": omega,
        "period_travel": period,
        "state_code": state,
    }
