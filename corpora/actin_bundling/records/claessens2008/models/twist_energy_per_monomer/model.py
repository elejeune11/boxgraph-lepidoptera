def run(inputs, settings):
    # d in nm -> m: 1e-9; tau in N*m^2; Delta_theta in rad -> result in N*m = J
    return {"E_twist": inputs["tau"] * (inputs["Delta_theta"] ** 2) / (inputs["d"] * 1e-9)}


def symbolic():
    import sympy as sp

    tau, Delta_theta, d = sp.symbols("tau Delta_theta d", positive=True)
    return {"E_twist": tau * Delta_theta**2 / (d * sp.Rational(1, 10**9))}
