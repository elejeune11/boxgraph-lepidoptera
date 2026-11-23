def run(inputs, settings):
    # Convert to SI (m, s, Pa)
    tau = inputs["tau"]
    D = inputs["D"] * 1e-12
    sigma_a = inputs["sigma_a"]
    xi = inputs["xi"]
    eta = inputs["eta"]
    q = inputs["q"] * 1e6
    L2 = eta / xi
    lam = -1.0 / tau - D * q * q + (sigma_a / xi) * (q * q) / (1.0 + L2 * q * q)
    return {"lam": lam}


def symbolic():
    import sympy as sp

    tau, D, sigma_a, xi, eta, q = sp.symbols("tau D sigma_a xi eta q", positive=True)
    D_si = D * sp.Rational(1, 10**12)
    q_si = q * 10**6
    L2 = eta / xi
    return {"lam": -1 / tau - D_si * q_si**2 + (sigma_a / xi) * q_si**2 / (1 + L2 * q_si**2)}
