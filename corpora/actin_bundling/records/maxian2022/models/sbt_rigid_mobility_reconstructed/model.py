import numpy as np
from numpy.polynomial.legendre import legvander


def run(inputs, settings):
    eps = float(inputs["epsilon"])
    delta = float(inputs["delta_rel"])
    N = int(round(settings.get("N", 16)))

    # Chebyshev Type-1 interior collocation nodes on [-1, 1], mapped to s in [0, 1]
    k = np.arange(N)
    u = -np.cos((2 * k + 1) * np.pi / (2 * N))
    V_leg = legvander(u, N - 1)
    rhs_w = np.zeros(N)
    rhs_w[0] = 2.0
    w = np.linalg.solve(V_leg.T, rhs_w)
    s = 0.5 * (u + 1.0)
    ws = 0.5 * w

    # Exact Legendre finite-part integral operator:
    # int_0^1 (P_n(s') - P_n(s)) / |s - s'| ds' = -2 * H_n * P_n(s)
    H = np.array([sum(1.0 / j for j in range(1, n + 1)) for n in range(N)])
    FP = V_leg @ np.diag(-2.0 * H) @ np.linalg.inv(V_leg)

    # C^1 ellipsoidal-cap endpoint regularization of c(s) = ln(4 s(1-s)/eps^2) over distance delta
    arg = np.where(
        s < delta,
        4.0 * (1.0 - s) * (delta**2) / ((eps**2) * (2.0 * delta - s)),
        np.where(
            s > 1.0 - delta,
            4.0 * s * (delta**2) / ((eps**2) * (2.0 * delta - (1.0 - s))),
            4.0 * s * (1.0 - s) / (eps**2),
        ),
    )
    c = np.log(arg)

    # Straight-rod perpendicular and parallel 1D mobility operators (Eq. 4, nondimensionalized by mu*L)
    M_perp = (np.diag(c + 1.0) + FP) / (8.0 * np.pi)
    M_par = (np.diag(2.0 * c - 2.0) + 2.0 * FP) / (8.0 * np.pi)

    # Rigid-body Schur complement N = (K_r^* M^-1 K_r)^dagger (Eq. 11)
    ones = np.ones(N)
    x_rot = s - 0.5
    alpha = float(1.0 / (ws @ np.linalg.solve(M_perp, ones)))
    alpha_plus_beta = float(1.0 / (ws @ np.linalg.solve(M_par, ones)))
    beta = alpha_plus_beta - alpha
    gamma = float(1.0 / ((ws * x_rot) @ np.linalg.solve(M_perp, x_rot)))
    return {"alpha": alpha, "beta": beta, "gamma": gamma}
