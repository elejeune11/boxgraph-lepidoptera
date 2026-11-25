def run(inputs, settings):
    return {"v_M": inputs["v_M_max"] * (1.0 + inputs["f"] / inputs["f_M_stall"])}

def symbolic():
    import sympy as sp
    v_M_max, f, f_M_stall = sp.symbols("v_M_max f f_M_stall", real=True)
    return {"v_M": v_M_max * (1 + f / f_M_stall)}
