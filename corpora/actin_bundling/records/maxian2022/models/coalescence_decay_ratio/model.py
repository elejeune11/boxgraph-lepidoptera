def run(inputs, settings):
    return {"ratio": inputs["tau_dec_nb"] / inputs["tau_dec_b"]}


def symbolic():
    import sympy as sp

    tau_dec_nb, tau_dec_b = sp.symbols("tau_dec_nb tau_dec_b", positive=True)
    return {"ratio": tau_dec_nb / tau_dec_b}
