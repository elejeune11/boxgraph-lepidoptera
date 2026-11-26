def run(inputs, settings):
    return {"P": inputs["zdm0"] / (inputs["gamma"] * inputs["D"])}


def symbolic():
    import sympy as sp

    zdm0, gamma, D = sp.symbols("zdm0 gamma D", positive=True)
    return {"P": zdm0 / (gamma * D)}
