def run(inputs, settings):
    return {"e_d": abs(inputs["x_bs"] - inputs["x_cg"]) / inputs["L"]}

def symbolic():
    import sympy as sp
    x_bs, x_cg, L = sp.symbols("x_bs x_cg L", real=True)
    return {"e_d": sp.Abs(x_bs - x_cg) / L}
