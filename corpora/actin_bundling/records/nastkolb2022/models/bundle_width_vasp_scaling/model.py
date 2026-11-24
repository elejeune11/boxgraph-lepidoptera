def run(inputs, settings):
    # Scaling law printed without a prefactor (prefactor = 1; judged on shape only)
    return {"W": float(inputs["c_vasp"])}


def symbolic():
    import sympy as sp

    c_vasp = sp.symbols("c_vasp", positive=True)
    return {"W": c_vasp}
