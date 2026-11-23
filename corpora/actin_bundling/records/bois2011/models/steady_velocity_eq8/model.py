def run(inputs, settings):
    return {"v": inputs["D"] * inputs["du_dx"]}


def symbolic():
    import sympy as sp

    D = sp.Symbol("D", positive=True)
    du_dx = sp.Symbol("du_dx", real=True)
    return {"v": D * du_dx}
