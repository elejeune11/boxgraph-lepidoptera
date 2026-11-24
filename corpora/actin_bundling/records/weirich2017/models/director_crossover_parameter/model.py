def run(inputs, settings):
    # gamma_A [nN/m = 1e-9 N/m] * V^(1/3) [um = 1e-6 m] = 1e-15 N = 1e-3 pN
    return {"chi": inputs["gamma_A"] * (inputs["V"] ** (1.0 / 3.0)) * 1e-3 / inputs["K"]}


def symbolic():
    import sympy as sp

    gamma_A, V, K = sp.symbols("gamma_A V K", positive=True)
    return {"chi": sp.Rational(1, 1000) * gamma_A * V ** sp.Rational(1, 3) / K}
