def run(inputs, settings):
    return {"U_f": inputs["U_f_stretch"] + inputs["U_f_bend"] + inputs["U_f_confine"]}


def symbolic():
    import sympy as sp

    U_f_stretch, U_f_bend, U_f_confine = sp.symbols("U_f_stretch U_f_bend U_f_confine", real=True)
    return {"U_f": U_f_stretch + U_f_bend + U_f_confine}
