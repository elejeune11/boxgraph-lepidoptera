def run(inputs, settings):
    return {"c_total": inputs["c_unlabeled"] + inputs["c_labeled"]}


def symbolic():
    import sympy as sp
    c_unlabeled, c_labeled = sp.symbols("c_unlabeled c_labeled", positive=True)
    return {"c_total": c_unlabeled + c_labeled}
