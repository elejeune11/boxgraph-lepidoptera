def run(inputs, settings):
    return {"nu": 1.0 + inputs["eta_v_over_eta"]}


def symbolic():
    import sympy as sp

    eta_v_over_eta = sp.symbols("eta_v_over_eta", positive=True)
    return {"nu": 1 + eta_v_over_eta}
