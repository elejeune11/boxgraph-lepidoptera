def run(inputs, settings):
    return {"R": inputs["l_coop"] / inputs["l_simple"]}


def symbolic():
    import sympy as sp

    l_coop, l_simple = sp.symbols("l_coop l_simple", positive=True)
    return {"R": l_coop / l_simple}
