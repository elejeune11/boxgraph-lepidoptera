import math
import numpy as np


def run(inputs, settings):
    kappa = float(inputs["kappa"])
    sigma = float(inputs["sigma"])
    kBT = float(inputs["kBT"])
    Lambda = float(inputs["Lambda"])
    H0 = float(inputs["H0"])
    ns = float(inputs["ns"])
    n0 = float(inputs["n0"])
    J = float(inputs["J"])
    xi = float(inputs["xi"])
    A = float(inputs["A"])
    L = float(inputs["L"])
    amp = float(inputs["n_pert_amp"])
    t_end = float(inputs["t_end"])

    nx = int(round(settings.get("nx", 64)))
    dt = float(settings.get("dt", 0.05))
    n_steps = max(1, int(round(t_end / dt)))
    dt = t_end / n_steps

    x = np.linspace(0.0, L, nx, endpoint=False)
    k = 2.0 * np.pi * np.fft.rfftfreq(nx, d=L / nx)
    k2 = k * k
    k3 = k2 * k
    k4 = k2 * k2
    dealias = np.abs(k) <= (2.0 / 3.0) * np.max(np.abs(k))

    D = Lambda * kBT
    c_disp = kappa * Lambda * (H0 ** 2) / ns
    c_curv = kappa * Lambda * H0 / ns
    c_agg = 2.0 * J * Lambda / ns

    # Linear stiff diagonal operators treated implicitly (IMEX-RK2 / Heun-Crank-Nicolson)
    Lh = -(kappa * k4 + sigma * k2) / xi
    Ln = -(D + c_disp * n0) * k2 - (c_agg * n0 / ns) * k4

    def rhs_nl(h_hat, n_hat):
        h_phys_x3 = np.fft.irfft(-1j * k3 * h_hat * dealias, n=nx)
        n_phys = np.fft.irfft(n_hat * dealias, n=nx)
        n_x = np.fft.irfft(1j * k * n_hat * dealias, n=nx)
        n_x3 = np.fft.irfft(-1j * k3 * n_hat * dealias, n=nx)

        # Residual nonlinear flux (after removing linear diagonal parts already in Ln):
        # J_tot - J_lin_diag = -c_disp*(n - n0)*n_x + c_agg*n*n_x + (c_agg/ns)*(n - n0)*n_x3 + c_curv*n*h_x3
        J_nl = (
            -c_disp * (n_phys - n0) * n_x
            + c_agg * n_phys * n_x
            + (c_agg / ns) * (n_phys - n0) * n_x3
            + c_curv * n_phys * h_phys_x3
        )
        nl_n = -1j * k * np.fft.rfft(J_nl) * dealias
        dn_hat = n_hat.copy()
        dn_hat[0] = 0.0
        nl_h = ((A - kappa * H0 * k2) / xi) * dn_hat
        return nl_h, nl_n

    # Deterministic multi-mode initial condition seeded by mode 4 (~q_max) + mode 2 + mode 5
    h = np.zeros(nx)
    n = n0 + amp * (
        np.cos(2.0 * np.pi * 4.0 * x / L)
        + 0.3 * np.cos(2.0 * np.pi * 2.0 * x / L + 0.5)
        + 0.2 * np.cos(2.0 * np.pi * 5.0 * x / L + 1.1)
    )
    h_hat = np.fft.rfft(h)
    n_hat = np.fft.rfft(n)

    # IMEX second-order Crank-Nicolson / Heun predictor-corrector
    denom_h = 1.0 - 0.5 * dt * Lh
    numer_h = 1.0 + 0.5 * dt * Lh
    denom_n = 1.0 - 0.5 * dt * Ln
    numer_n = 1.0 + 0.5 * dt * Ln

    step_t1 = max(1, int(round(0.3 * n_steps)))
    step_t2 = max(step_t1 + 1, int(round(0.6 * n_steps)))
    amp_m4_t1 = None
    amp_m4_t2 = None

    for step in range(1, n_steps + 1):
        nl_h1, nl_n1 = rhs_nl(h_hat, n_hat)
        # Predictor
        h_pred = (numer_h * h_hat + dt * nl_h1) / denom_h
        n_pred = (numer_n * n_hat + dt * nl_n1) / denom_n
        # Corrector
        nl_h2, nl_n2 = rhs_nl(h_pred, n_pred)
        h_hat = (numer_h * h_hat + 0.5 * dt * (nl_h1 + nl_h2)) / denom_h
        n_hat = (numer_n * n_hat + 0.5 * dt * (nl_n1 + nl_n2)) / denom_n
        if step == step_t1:
            amp_m4_t1 = abs(n_hat[4]) / nx
        elif step == step_t2:
            amp_m4_t2 = abs(n_hat[4]) / nx

    h_final = np.fft.irfft(h_hat, n=nx)
    n_final = np.fft.irfft(n_hat, n=nx)
    mode4_rate = math.log(amp_m4_t2 / amp_m4_t1) / ((step_t2 - step_t1) * dt)

    def _subgrid_ext(arr, find_max=True):
        idx = int(np.argmax(arr) if find_max else np.argmin(arr))
        y0 = float(arr[(idx - 1) % nx])
        y1 = float(arr[idx])
        y2 = float(arr[(idx + 1) % nx])
        denom = y0 - 2.0 * y1 + y2
        if abs(denom) > 1e-30:
            delta = 0.5 * (y0 - y2) / denom
            if abs(delta) <= 1.0:
                return y1 - 0.25 * (y0 - y2) * delta
        return y1

    return {
        "n_max": float(_subgrid_ext(n_final, True)),
        "n_min": float(_subgrid_ext(n_final, False)),
        "h_max": float(_subgrid_ext(h_final, True)),
        "mode4_growth_rate": float(mode4_rate),
    }
