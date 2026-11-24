def run(inputs, settings):
    # (L_target - L0) in um -> * 1000 nm; divided by T_grow in s -> nm/s
    return {"k_grow": 1000.0 * (inputs["L_target"] - inputs["L0"]) / inputs["T_grow"]}


def symbolic():
    import sympy as sp

    L_target, L0, T_grow = sp.symbols("L_target L0 T_grow", positive=True)
    return {"k_grow": 1000 * (L_target - L0) / T_grow}
