def run(inputs, settings):
    # Each fascin bridges 2 filaments, so z binding sites around a filament use z/2 fascins per filament
    return {"u": 0.5 * inputs["z"] * inputs["r_sat"]}


def symbolic():
    import sympy as sp

    z, r_sat = sp.symbols("z r_sat", positive=True)
    return {"u": sp.Rational(1, 2) * z * r_sat}
