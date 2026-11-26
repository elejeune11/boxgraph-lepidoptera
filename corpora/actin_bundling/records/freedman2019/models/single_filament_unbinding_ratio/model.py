def run(inputs, settings):
    return {"off1_ratio": inputs["k_off_1_short"] / inputs["k_off_1_long"]}


def symbolic():
    import sympy as sp

    k_off_1_short, k_off_1_long = sp.symbols("k_off_1_short k_off_1_long", positive=True)
    return {"off1_ratio": k_off_1_short / k_off_1_long}
