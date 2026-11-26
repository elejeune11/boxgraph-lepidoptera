"""Cross-link finite-thickness geometric factor chi^2 = A_f * (d_f + t)^2 / (12 * I_f)."""


def run(inputs, settings):
    return {"chi2": inputs["A_f"] * ((inputs["d_f"] + inputs["t"]) ** 2) / (12.0 * inputs["I_f"])}


def symbolic():
    import sympy as sp

    A_f, d_f, t, I_f = sp.symbols("A_f d_f t I_f", positive=True)
    return {"chi2": A_f * (d_f + t) ** 2 / (12 * I_f)}
