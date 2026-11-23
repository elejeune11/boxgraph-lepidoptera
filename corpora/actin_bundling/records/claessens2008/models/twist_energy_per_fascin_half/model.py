K_B = 1.380649e-23  # J/K


def run(inputs, settings):
    E_mon = 0.5 * inputs["tau"] * (inputs["Delta_theta"] ** 2) / (inputs["d"] * 1e-9)
    kBT = K_B * inputs["T"]
    return {"E_twist_fascin_kBT": E_mon / (inputs["R_star"] * kBT)}


def symbolic():
    import sympy as sp

    tau, Delta_theta, d, R_star, T = sp.symbols("tau Delta_theta d R_star T", positive=True)
    return {
        "E_twist_fascin_kBT": sp.Rational(1, 2)
        * (tau * Delta_theta**2 / (d * sp.Rational(1, 10**9)))
        / (R_star * sp.Float(K_B) * T)
    }
