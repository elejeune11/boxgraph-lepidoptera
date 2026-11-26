def run(inputs, settings):
    return {"v_inv": inputs["tau"] / inputs["d"]}


def symbolic():
    import sympy as sp

    tau, d = sp.symbols("tau d", positive=True)
    return {"v_inv": tau / d}
