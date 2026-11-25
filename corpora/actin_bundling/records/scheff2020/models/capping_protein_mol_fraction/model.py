def run(inputs, settings):
    # c_cp [nM] / (c_actin [uM] * 1000 [nM/uM]) * 100 [%] = c_cp / (10 * c_actin)
    return {"f_cp": inputs["c_cp"] / (10.0 * inputs["c_actin"])}


def symbolic():
    import sympy as sp

    c_actin, c_cp = sp.symbols("c_actin c_cp", positive=True)
    return {"f_cp": c_cp / (10 * c_actin)}
