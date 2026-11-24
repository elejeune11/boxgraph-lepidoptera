def run(inputs, settings):
    return {"V1_over_V2": inputs["Q_1"] / inputs["Q_2"]}


def symbolic():
    import sympy as sp

    Q_1, Q_2 = sp.symbols("Q_1 Q_2", positive=True)
    return {"V1_over_V2": Q_1 / Q_2}
