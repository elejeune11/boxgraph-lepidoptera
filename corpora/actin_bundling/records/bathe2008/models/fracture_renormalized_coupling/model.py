"""Fracture-renormalized coupling parameter: alpha' = alpha * (phi* / phi)^2."""


def run(inputs, settings):
    return {"alpha_prime": inputs["alpha"] * ((inputs["phi_star"] / inputs["phi"]) ** 2)}


def symbolic():
    import sympy as sp

    alpha, phi_star, phi = sp.symbols("alpha phi_star phi", positive=True)
    return {"alpha_prime": alpha * (phi_star / phi) ** 2}
