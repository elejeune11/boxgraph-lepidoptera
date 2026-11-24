import math
from scipy.optimize import root


def run(inputs, settings):
    beta = float(inputs["beta"])
    E1 = float(inputs["E1"])
    E3 = float(inputs["E3"])
    tol = float(settings.get("tol", 1e-12))

    # Effective bulk free energy after minimizing Eq. (10)'s bulk term Q_kl Q_kl = max(0, (1+phi)/beta):
    # f_eff(phi) = phi^4/4 - phi^2/2 - (1+phi)^2 / (4*beta) for phi > -1, else phi^4/4 - phi^2/2.
    def f_eff(p):
        base = 0.25 * p**4 - 0.5 * p**2
        return base - (1.0 + p) ** 2 / (4.0 * beta) if p > -1.0 else base

    def mu_eff(p):
        base = p**3 - p
        return base - (1.0 + p) / (2.0 * beta) if p > -1.0 else base

    def residual(x):
        p1, p2 = float(x[0]), float(x[1])
        m1, m2 = mu_eff(p1), mu_eff(p2)
        g1, g2 = f_eff(p1) - m1 * p1, f_eff(p2) - m2 * p2
        return [m1 - m2, g1 - g2]

    sol = root(residual, [-1.2, 1.25], tol=tol)
    if not sol.success:
        raise RuntimeError(f"Common-tangent solver failed: {sol.message}")
    p_lo, p_ld = sorted([float(sol.x[0]), float(sol.x[1])])
    mu_coex = mu_eff(p_lo)
    Q_eq = math.sqrt(max(0.0, (1.0 + p_ld) / (2.0 * beta)))
    eff_diff_perp = E1 - E3 * Q_eq
    return {
        "phi_Lo_eq": p_lo,
        "phi_Ld_eq": p_ld,
        "mu_coex": mu_coex,
        "Q_eq_Ld": Q_eq,
        "eff_diff_perp": eff_diff_perp,
    }
