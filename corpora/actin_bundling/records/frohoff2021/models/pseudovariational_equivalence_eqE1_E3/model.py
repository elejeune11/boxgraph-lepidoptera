import numpy as np


def run(inputs, settings):
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    Q = inputs["Q"]
    alpha = inputs["alpha"]
    ell = inputs["ell"]
    amp1 = inputs["amp1"]
    amp2 = inputs["amp2"]
    Nx = int(settings.get("n_grid", 64))

    x = np.linspace(0.0, 1.0, Nx, endpoint=False)
    # Generic smooth test fields on periodic [-1/2, 1/2] with relative phase shift
    p1 = amp1 * (np.cos(2.0 * np.pi * x) + 0.25 * np.sin(4.0 * np.pi * x))
    p2 = amp2 * (np.sin(2.0 * np.pi * x) + 0.2 * np.cos(4.0 * np.pi * x))

    k = 2.0 * np.pi * np.fft.rfftfreq(Nx, d=1.0 / Nx)
    ik = 1j * k
    k2 = k * k

    def d1(u):
        return np.fft.irfft(ik * np.fft.rfft(u), Nx)

    def d2(u):
        return np.fft.irfft(-k2 * np.fft.rfft(u), Nx)

    # Eq. (4) RHS at rho = 0
    mu1_eq4 = -(1.0 / ell**2) * d2(p1) + a * p1 + p1**3 - alpha * p2
    mu2_eq4 = -(kappa / ell**2) * d2(p2) + (a + a_Delta) * p2 + p2**3 + alpha * p1
    rhs1_eq4 = (1.0 / ell**2) * d2(mu1_eq4)
    rhs2_eq4 = (Q / ell**2) * d2(mu2_eq4)

    # Eqs. (E1)-(E3): Q1 = -1, Q2 = Q
    dF_dp1 = (1.0 / ell**2) * d2(p1) - (a * p1 + p1**3) + alpha * p2
    dF_dp2 = -(kappa / ell**2) * d2(p2) + ((a + a_Delta) * p2 + p2**3) + alpha * p1
    Q1 = -1.0
    Q2 = Q
    rhs1_eqE = (Q1 / ell**2) * d2(dF_dp1)
    rhs2_eqE = (Q2 / ell**2) * d2(dF_dp2)

    rhs_diff = float(max(np.max(np.abs(rhs1_eq4 - rhs1_eqE)), np.max(np.abs(rhs2_eq4 - rhs2_eqE))))
    # d(F_tilde)/dt = -int [ (Q1/ell^2)*(d_x dF_dp1)^2 + (Q2/ell^2)*(d_x dF_dp2)^2 ] dx
    dF_dt = float(-np.mean((Q1 / ell**2) * d1(dF_dp1) ** 2 + (Q2 / ell**2) * d1(dF_dp2) ** 2))

    return {"rhs_max_diff": rhs_diff, "dF_tilde_dt": dF_dt}
