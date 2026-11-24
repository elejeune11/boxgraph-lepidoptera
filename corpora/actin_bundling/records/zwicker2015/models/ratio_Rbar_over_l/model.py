def run(inputs, settings):
    return {"Rbar_over_l": inputs["R_bar"] / inputs["l"]}


def symbolic():
    import sympy as sp

    R_bar, l = sp.symbols("R_bar l", positive=True)
    return {"Rbar_over_l": R_bar / l}
