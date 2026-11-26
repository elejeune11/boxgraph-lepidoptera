def run(inputs, settings):
    f_i = inputs["F_xl"] * (inputs["d_next"] / inputs["l_link"])
    f_next = inputs["F_xl"] - f_i
    return {"F_i": f_i, "F_i_plus_1": f_next}


def symbolic():
    import sympy as sp

    F_xl, d_next, l_link = sp.symbols("F_xl d_next l_link", positive=True)
    f_i = F_xl * d_next / l_link
    return {"F_i": f_i, "F_i_plus_1": F_xl - f_i}
