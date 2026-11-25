def run(inputs, settings):
    Delta_g = float(inputs["Delta_g"])
    Delta_psi = float(inputs["Delta_psi"])
    sigma = float(inputs["sigma"])
    d = float(inputs["d"])
    R = float(inputs["R"])
    return {"mu_I": Delta_g / Delta_psi + (d - 1.0) * sigma / (Delta_psi * R)}


def symbolic():
    import sympy as sp

    Delta_g, Delta_psi, sigma, d, R = sp.symbols("Delta_g Delta_psi sigma d R", real=True)
    return {"mu_I": Delta_g / Delta_psi + (d - 1) * sigma / (Delta_psi * R)}
