def run(inputs, settings):
    return {"E_twist": 0.5 * inputs["tau"] * (inputs["Delta_theta"] ** 2) / (inputs["d"] * 1e-9)}


def symbolic():
    import sympy as sp

    tau, Delta_theta, d = sp.symbols("tau Delta_theta d", positive=True)
    return {"E_twist": sp.Rational(1, 2) * tau * Delta_theta**2 / (d * sp.Rational(1, 10**9))}
