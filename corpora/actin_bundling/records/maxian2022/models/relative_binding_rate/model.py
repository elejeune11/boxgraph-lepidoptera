def run(inputs, settings):
    return {"kon_rel": inputs["kon"] / inputs["kon0"]}


def symbolic():
    import sympy as sp

    kon, kon0 = sp.symbols("kon kon0", positive=True)
    return {"kon_rel": kon / kon0}
