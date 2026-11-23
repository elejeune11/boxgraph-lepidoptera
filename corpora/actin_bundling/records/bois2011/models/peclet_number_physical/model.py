def run(inputs, settings):
    return {"Pe": inputs["zdm0"] / (inputs["D"] * inputs["gamma"])}


def symbolic():
    import sympy as sp

    zdm0, D, gamma = sp.symbols("zdm0 D gamma", positive=True)
    return {"Pe": zdm0 / (D * gamma)}
