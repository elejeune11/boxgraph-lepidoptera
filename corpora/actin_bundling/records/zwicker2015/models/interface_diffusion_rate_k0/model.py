def run(inputs, settings):
    return {"k_0": inputs["D_B"] / (inputs["w"] ** 2)}


def symbolic():
    import sympy as sp

    D_B, w = sp.symbols("D_B w", positive=True)
    return {"k_0": D_B / w**2}
