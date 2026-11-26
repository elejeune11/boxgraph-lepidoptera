def run(inputs, settings):
    return {"c_loc": inputs["theta_18"] / inputs["V_P"]}


def symbolic():
    import sympy as sp

    theta_18, V_P = sp.symbols("theta_18 V_P", positive=True)
    return {"c_loc": theta_18 / V_P}
