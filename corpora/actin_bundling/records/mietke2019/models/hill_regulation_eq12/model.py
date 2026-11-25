def run(inputs, settings):
    c = inputs["c"]
    c0 = inputs["c0"]
    m = inputs["m"]
    cm = c ** m
    c0m = c0 ** m
    return {
        "f": cm / (c0m + cm),
        "c0_df_dc0": m / 4.0,
    }


def symbolic():
    import sympy as sp

    c, c0, m = sp.symbols("c c0 m", positive=True)
    f_expr = c**m / (c0**m + c**m)
    df_c0 = sp.simplify(c0 * sp.diff(f_expr, c).subs(c, c0))
    return {"f": f_expr, "c0_df_dc0": df_c0}
