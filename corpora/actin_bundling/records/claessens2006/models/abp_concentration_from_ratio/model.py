def run(inputs, settings):
    return {"c_ABP": inputs["r_ABP_a"] * inputs["c_a"]}


def symbolic():
    import sympy as sp

    r_ABP_a, c_a = sp.symbols("r_ABP_a c_a", positive=True)
    return {"c_ABP": r_ABP_a * c_a}
