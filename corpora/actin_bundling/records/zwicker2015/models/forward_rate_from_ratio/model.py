def run(inputs, settings):
    return {"k_f": inputs["kf_over_kb"] * inputs["k_b_in"]}


def symbolic():
    import sympy as sp

    k_b_in, kf_over_kb = sp.symbols("k_b_in kf_over_kb", positive=True)
    return {"k_f": kf_over_kb * k_b_in}
