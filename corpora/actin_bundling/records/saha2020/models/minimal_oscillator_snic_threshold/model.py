from scipy.optimize import brentq, minimize_scalar


def run(inputs, settings):
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    chi = float(inputs["chi"])
    tol = float(settings.get("tol", 1e-10))

    # For alpha > |chi|, nontrivial fixed points satisfy:
    #   x2 = -4 * x1 * (x1^2 - c1^2) / (chi + alpha)
    #   4 * g(u) * (u * g(u)^2 - c2^2) + (chi - alpha) = 0, with u = x1^2 in (0, c1^2)
    # where g(u) = -4 * (u - c1^2) / (chi + alpha) > 0 when chi < 0 and alpha > |chi|.
    # Handle general sign of chi via symmetry alpha -> alpha, chi -> -|chi| (sign of chi flips x2 -> -x2).
    chi_neg = -abs(chi)

    def max_F(alpha):
        denom = chi_neg + alpha
        if denom <= 1e-14:
            return 1e6

        def neg_F(u):
            g = -4.0 * (u - c1 * c1) / denom
            return -(4.0 * g * (u * g * g - c2 * c2) + (chi_neg - alpha))

        res = minimize_scalar(neg_F, bounds=(1e-12 * c1 * c1, (1.0 - 1e-12) * c1 * c1), method="bounded")
        return -float(res.fun)

    lo = abs(chi) + 1e-6
    hi = abs(chi) + 0.05
    while max_F(hi) > 0.0:
        hi *= 2.0
    alpha_sn = brentq(max_F, lo, hi, xtol=tol, rtol=tol)
    return {"alpha_sn": float(alpha_sn)}
