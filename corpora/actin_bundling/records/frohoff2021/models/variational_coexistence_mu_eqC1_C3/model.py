import math
from scipy.optimize import root


def run(inputs, settings):
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    rho = inputs["rho"]
    mu1_target = inputs["mu1"]
    branch = str(settings.get("branch", "I-IV"))
    tol = float(settings.get("tol", 1e-10))

    def chem_and_omega(p1, p2):
        mu1 = a * p1 + p1**3 - rho * p2
        mu2 = (a + a_Delta) * p2 + p2**3 - rho * p1
        f1 = 0.5 * a * p1**2 + 0.25 * p1**4
        f2 = 0.5 * (a + a_Delta) * p2**2 + 0.25 * p2**4
        w = f1 + f2 - rho * p1 * p2 - mu1 * p1 - mu2 * p2
        return mu1, mu2, w

    def res(x):
        p1A, p2A, p1B, p2B = x
        mu1A, mu2A, wA = chem_and_omega(p1A, p2A)
        mu1B, mu2B, wB = chem_and_omega(p1B, p2B)
        return [mu1A - mu1_target, mu1B - mu1_target, mu2A - mu2B, wA - wB]

    p1_g = math.sqrt(max(0.5, -a))
    p2_g = math.sqrt(max(0.5, -(a + a_Delta)))
    if branch == "I-IV":
        # High-phi1 phases I (+,+) and IV (+,-)
        x0 = [max(p1_g, 1.5), p2_g, max(p1_g, 1.4), -p2_g]
    elif branch == "I-III":
        # Diagonal coexistence I (+,+) and III (-,-)
        x0 = [p1_g + 0.3, p2_g + 0.3, -p1_g - 0.3, -p2_g - 0.3]
    elif branch == "II-III":
        # Low-phi1 phases II (-,+) and III (-,-)
        x0 = [-max(p1_g, 1.4), p2_g, -max(p1_g, 1.5), -p2_g]
    else:
        # I-II: (+,+) and (-,+)
        x0 = [p1_g, max(p2_g, 1.5), -p1_g, max(p2_g, 1.4)]

    sol = root(res, x0, method="hybr", tol=tol)
    if not sol.success and max(abs(v) for v in sol.fun) > max(1e-8, 100.0 * tol):
        return {"mu2": None, "phi1_A": None, "phi2_A": None}
    p1A, p2A, p1B, p2B = [float(v) for v in sol.x]
    _, mu2, _ = chem_and_omega(p1A, p2A)
    return {"mu2": float(mu2), "phi1_A": p1A, "phi2_A": p2A}
