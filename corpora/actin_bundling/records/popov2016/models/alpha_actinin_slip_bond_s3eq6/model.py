import math

def run(inputs, settings):
    return {"k_unbind": inputs["k_unbind_0"] * math.exp(inputs["F_ext"] / inputs["F_unbind"])}

def symbolic():
    import sympy as sp
    k_unbind_0, F_ext, F_unbind = sp.symbols("k_unbind_0 F_ext F_unbind", positive=True)
    return {"k_unbind": k_unbind_0 * sp.exp(F_ext / F_unbind)}
