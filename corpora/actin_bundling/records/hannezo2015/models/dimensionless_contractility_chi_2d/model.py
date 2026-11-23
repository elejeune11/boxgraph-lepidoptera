def run(inputs, settings):
    # sigma_a [Pa*m] / (xi [Pa*s/m] * D [um^2/s -> 1e-12 m^2/s]) -> dimensionless
    return {"chi": inputs["sigma_a"] / (inputs["xi"] * inputs["D"] * 1e-12)}


def symbolic():
    import sympy as sp

    sigma_a, xi, D = sp.symbols("sigma_a xi D", positive=True)
    return {"chi": sigma_a / (xi * D * sp.Rational(1, 10**12))}
