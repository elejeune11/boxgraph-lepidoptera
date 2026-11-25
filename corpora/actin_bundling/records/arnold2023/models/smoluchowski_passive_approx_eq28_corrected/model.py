import math

K_B = 1.380649e-23  # J / K


def run(inputs, settings):
    T = float(inputs["T"])
    eta = float(inputs["eta"])
    phi_A = float(inputs["phi_A"])
    L = float(inputs["L"])
    a0 = float(inputs["a_0"])
    eps = float(inputs["epsilon"])
    t = float(inputs["t"])

    denom = 3.0 * eps + 1.0
    k = (9.0 * K_B * T * phi_A / (16.0 * eta * denom)) * 1.0e18  # um^3 / s
    b = (3.0 * math.log(L / (2.0 * a0)) + 1.0) / denom

    arg = k * t + b * (a0 ** 3)
    a = arg ** (1.0 / 3.0)
    da_dt = k / (3.0 * a * a)
    alpha_eff = (t / a) * da_dt
    return {
        "k": float(k),
        "b": float(b),
        "a": float(a),
        "da_dt": float(da_dt),
        "alpha_eff": float(alpha_eff),
    }


def symbolic():
    import sympy as sp
    T, eta, phi_A, L, a_0, epsilon, t = sp.symbols("T eta phi_A L a_0 epsilon t", positive=True)
    denom = 3 * epsilon + 1
    k = (9 * K_B * T * phi_A / (16 * eta * denom)) * 1e18
    b = (3 * sp.log(L / (2 * a_0)) + 1) / denom
    a = (k * t + b * a_0**3) ** sp.Rational(1, 3)
    da_dt = k / (3 * a**2)
    alpha_eff = (t / a) * da_dt
    return {"k": k, "b": b, "a": a, "da_dt": da_dt, "alpha_eff": alpha_eff}
