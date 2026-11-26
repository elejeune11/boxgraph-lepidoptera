def run(inputs, settings):
    a_0 = inputs["a_0"]
    a = inputs["a"]
    inv_lam0 = inputs["inv_lambda_0"]
    q = inputs["q"]
    f_pack = inputs["f_pack"]
    # lambda_0 - q * alpha_b / v_eff = 0 => v_eff = q * pi * (a^2 - a_0^2) * inv_lambda_0
    # v_eff = (4*pi/3) * delta_eff^3 / f_pack => delta_eff^3 = (3/4) * f_pack * q * (a^2 - a_0^2) * inv_lambda_0
    return {"delta_eff_min": (0.75 * f_pack * q * (a * a - a_0 * a_0) * inv_lam0) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    a_0, a, inv_lambda_0, q, f_pack = sp.symbols("a_0 a inv_lambda_0 q f_pack", positive=True)
    return {"delta_eff_min": (sp.Rational(3, 4) * f_pack * q * (a**2 - a_0**2) * inv_lambda_0) ** sp.Rational(1, 3)}
