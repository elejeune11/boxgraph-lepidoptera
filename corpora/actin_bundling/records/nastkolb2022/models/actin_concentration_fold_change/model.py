def run(inputs, settings):
    return {"R_actin": inputs["c_actin_high"] / inputs["c_actin_std"]}


def symbolic():
    import sympy as sp

    c_actin_high, c_actin_std = sp.symbols("c_actin_high c_actin_std", positive=True)
    return {"R_actin": c_actin_high / c_actin_std}
