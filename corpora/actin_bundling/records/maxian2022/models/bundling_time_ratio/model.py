def run(inputs, settings):
    return {"ratio": inputs["tau_c_nb"] / inputs["tau_c_b"]}


def symbolic():
    import sympy as sp

    tau_c_nb, tau_c_b = sp.symbols("tau_c_nb tau_c_b", positive=True)
    return {"ratio": tau_c_nb / tau_c_b}
