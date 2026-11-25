def run(inputs, settings):
    return {"K_d": float(inputs["k_unbind"]) / float(inputs["k_bind_2"])}


def symbolic():
    import sympy as sp

    k_unbind, k_bind_2 = sp.symbols("k_unbind k_bind_2", positive=True)
    return {"K_d": k_unbind / k_bind_2}
