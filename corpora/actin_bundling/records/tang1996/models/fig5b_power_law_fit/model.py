def run(inputs, settings):
    return {"c_onset": inputs["A"] * (inputs["c_KCl"] ** inputs["alpha"])}


def symbolic():
    import sympy as sp

    A, c_KCl, alpha = sp.symbols("A c_KCl alpha", positive=True)
    return {"c_onset": A * (c_KCl ** alpha)}
