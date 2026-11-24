import math
from scipy.optimize import root


def run(inputs, settings):
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    rho = inputs["rho"]
    phi1_A = inputs["phi1_A"]
    tol = float(settings.get("tol", 1e-12))

    def chem_and_omega(p1, p2):
        mu1 = a * p1 + p1**3 - rho * p2
        mu2 = (a + a_Delta) * p2 + p2**3 - rho * p1
        f1 = 0.5 * a * p1**2 + 0.25 * p1**4
        f2 = 0.5 * (a + a_Delta) * p2**2 + 0.25 * p2**4
        w = f1 + f2 - rho * p1 * p2 - mu1 * p1 - mu2 * p2
        return mu1, mu2, w

    def res(x):
        p2A, p1B, p2B = x
        mu1A, mu2A, wA = chem_and_omega(phi1_A, p2A)
        mu1B, mu2B, wB = chem_and_omega(p1B, p2B)
        return [mu1A - mu1B, mu2A - mu2B, wA - wB]

    p2_guess = math.sqrt(max(0.1, -(a + a_Delta)))
    x0 = [p2_guess, phi1_A - 0.1 * rho, -p2_guess]
    sol = root(res, x0, method="hybr", tol=tol)
    if not sol.success and max(abs(v) for v in sol.fun) > max(1e-9, 10.0 * tol):
        return {
            "phi2_A": None,
            "phi1_B": None,
            "phi2_B": None,
            "mu1": None,
            "mu2": None,
            "omega": None,
        }
    p2A, p1B, p2B = [float(v) for v in sol.x]
    mu1, mu2, w = chem_and_omega(phi1_A, p2A)
    return {
        "phi2_A": p2A,
        "phi1_B": p1B,
        "phi2_B": p2B,
        "mu1": float(mu1),
        "mu2": float(mu2),
        "omega": float(w),
    }
