def run(inputs, settings):
    return {"k_b": inputs["kb_over_kf"] * inputs["k_f_in"]}


def symbolic():
    import sympy as sp

    k_f_in, kb_over_kf = sp.symbols("k_f_in kb_over_kf", positive=True)
    return {"k_b": kb_over_kf * k_f_in}
