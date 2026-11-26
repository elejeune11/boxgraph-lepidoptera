def run(inputs, settings):
    return {"ratio": inputs["Pi_coac"] / inputs["Pi_star"]}


def symbolic():
    import sympy as sp

    Pi_coac, Pi_star = sp.symbols("Pi_coac Pi_star", positive=True)
    return {"ratio": Pi_coac / Pi_star}
