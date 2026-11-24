def run(inputs, settings):
    return {"D_B": inputs["m_B"] * inputs["b"]}


def symbolic():
    import sympy as sp

    m_B, b = sp.symbols("m_B b", positive=True)
    return {"D_B": m_B * b}
