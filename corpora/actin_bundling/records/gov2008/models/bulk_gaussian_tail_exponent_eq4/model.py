import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    k_b_prime = inputs["k_b_prime"]
    rho = inputs["rho"]
    E_0 = inputs["E_0"]

    def phi(r):
        x = r / a
        return -math.pi * x * eps_b * (x - 2.0) / (x * x) + 2.0 * math.pi * k_b_prime * (a * a) * (r * r) * (rho * rho) * math.log(x) - E_0

    h = 1e-6
    p_lo, p_hi = phi(R * (1.0 - h)), phi(R * (1.0 + h))
    if p_lo <= 0.0 or p_hi <= 0.0:
        return {"tail_exponent": None}
    q = (math.log(p_hi) - math.log(p_lo)) / math.log((1.0 + h) / (1.0 - h))
    return {"tail_exponent": q}


def symbolic():
    import sympy as sp

    R, a, eps_b, k_b_prime, rho, E_0 = sp.symbols("R a eps_b k_b_prime rho E_0", positive=True)
    x = R / a
    phi = -sp.pi * x * eps_b * (x - 2) / x**2 + 2 * sp.pi * k_b_prime * a**2 * R**2 * rho**2 * sp.log(x) - E_0
    return {"tail_exponent": sp.simplify(R * sp.diff(phi, R) / phi)}
