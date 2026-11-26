def run(inputs, settings):
    # Scaling law printed without a prefactor (prefactor: unknown -> prefactor = 1).
    return {"msd": inputs["t"] ** inputs["nu"]}

def symbolic():
    import sympy as sp
    t, nu = sp.symbols("t nu", positive=True)
    return {"msd": t**nu}
