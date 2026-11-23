def run(inputs, settings):
    # sigma_a [Pa] / (D [um^2/s -> 1e-12 m^2/s] * chi_c) -> Pa*s/m^2
    return {"xi": inputs["sigma_a"] / (inputs["D"] * 1e-12 * inputs["chi_c"])}


def symbolic():
    import sympy as sp

    sigma_a, D, chi_c = sp.symbols("sigma_a D chi_c", positive=True)
    return {"xi": sigma_a / (D * sp.Rational(1, 10**12) * chi_c)}
