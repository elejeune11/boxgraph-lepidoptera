import numpy as np
from numpy.polynomial.legendre import legvander


def run(inputs, settings):
    eps = float(inputs["epsilon"])
    delta = float(inputs["delta_rel"])
    Kc = float(inputs["Kc"])
    mu = float(inputs["mu"])
    L = float(inputs["L"])
    N = int(round(settings.get("N", 50)))

    # Chebyshev Type-1 collocation nodes and polynomial quadrature weights on s in [0, 1]
    k = np.arange(N)
    u = -np.cos((2 * k + 1) * np.pi / (2 * N))
    V_leg = legvander(u, N - 1)
    rhs_w = np.zeros(N)
    rhs_w[0] = 2.0
    w = np.linalg.solve(V_leg.T, rhs_w)
    s = 0.5 * (u + 1.0)
    ws = 0.5 * w

    # C^1 ellipsoidal-cap endpoint regularization over delta
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

    # Local drag mobility (Eq. 4 without integral term, nondimensionalized by mu*L)
    M_perp = (c + 1.0) / (8.0 * np.pi)
    M_par = (2.0 * c - 2.0) / (8.0 * np.pi)

    alpha = float(1.0 / np.sum(ws / M_perp))
    alpha_plus_beta = float(1.0 / np.sum(ws / M_par))
    beta = alpha_plus_beta - alpha
    gamma = float(1.0 / np.sum(ws * ((s - 0.5) ** 2) / M_perp))

    # Relative perpendicular endpoint relaxation rate:
    # lambda_perp = 2 * (Kc / (mu * L)) * (alpha + gamma / 4)
    # 1 Pa*s = 1 pN*s/um^2 -> lambda_perp in 1/s
    lam_perp = 2.0 * (Kc / (mu * L)) * (alpha + 0.25 * gamma)
    dt_max = 2.0 / lam_perp
    return {"alpha": alpha, "beta": beta, "gamma": gamma, "dt_max": float(dt_max)}
