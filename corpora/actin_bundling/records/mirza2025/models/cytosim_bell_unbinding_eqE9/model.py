import math

def run(inputs, settings):
    return {"r_u": inputs["r_u0"] * math.exp(inputs["f"] / inputs["f_u"])}

def symbolic():
    import sympy as sp
    r_u0, f, f_u = sp.symbols("r_u0 f f_u", real=True)
    return {"r_u": r_u0 * sp.exp(f / f_u)}
