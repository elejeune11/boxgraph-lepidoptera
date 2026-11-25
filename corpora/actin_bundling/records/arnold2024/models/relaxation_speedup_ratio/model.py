def run(inputs, settings):
    return {"speedup": inputs["tau_Ld"] / inputs["tau_LoLd"]}


def symbolic():
    import sympy as sp

    tau_Ld, tau_LoLd = sp.symbols("tau_Ld tau_LoLd", positive=True)
    return {"speedup": tau_Ld / tau_LoLd}
