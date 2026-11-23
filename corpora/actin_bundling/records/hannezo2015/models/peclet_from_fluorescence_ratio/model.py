def run(inputs, settings):
    return {"chi_plateau": inputs["chi_c"] * (inputs["F_plateau"] / inputs["F_c"])}


def symbolic():
    import sympy as sp

    chi_c, F_c, F_plateau = sp.symbols("chi_c F_c F_plateau", positive=True)
    return {"chi_plateau": chi_c * F_plateau / F_c}
