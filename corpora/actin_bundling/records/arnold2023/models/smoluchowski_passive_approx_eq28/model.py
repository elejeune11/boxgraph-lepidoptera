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

    # Printed Eq. (28): denominator 3*ln(1/2) + 3*eps + 1 is NEGATIVE (~ -1.079 + 3*eps) for eps << 1!
    denom = 3.0 * math.log(0.5) + 3.0 * eps + 1.0
    k = (9.0 * K_B * T * phi_A / (16.0 * eta * denom)) * 1.0e18  # um^3 / s
    # Printed Eq. (28): b = 3*ln(L/(2*a0)) + 1 (not divided by denom)
    b = 3.0 * math.log(L / (2.0 * a0)) + 1.0

    arg = k * t + b * (a0 ** 3)
    a = math.copysign(abs(arg) ** (1.0 / 3.0), arg)
    da_dt = k / (3.0 * a * a) if abs(a) > 0.0 else 0.0
    # Signed local exponent (t / |a|) * (da/dt): equals exact d(ln a)/d(ln t) when a > 0, negative when k < 0
    alpha_eff = (t / max(abs(a), 1.0e-30)) * da_dt
    return {
        "k": float(k),
        "b": float(b),
        "a": float(a),
        "da_dt": float(da_dt),
        "alpha_eff": float(alpha_eff),
    }
