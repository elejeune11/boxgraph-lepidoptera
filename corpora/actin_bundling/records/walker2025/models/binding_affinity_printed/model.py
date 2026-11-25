def run(inputs, settings):
    return {"K_d_ratio": float(inputs["k_unbind"]) / float(inputs["k_bind"])}


def symbolic():
    import sympy as sp

    k_unbind, k_bind = sp.symbols("k_unbind k_bind", positive=True)
    return {"K_d_ratio": k_unbind / k_bind}
