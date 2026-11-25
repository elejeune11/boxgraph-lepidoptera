def run(inputs, settings):
    f = inputs["f_gamma"]
    return {
        "k_bar_d": inputs["k_bar_d_ref"] / f,
        "c_bar_0": inputs["c_bar_0_ref"] / f,
        "b_bar": inputs["b_bar_ref"] / f,
        "lambda_o_bar": inputs["lambda_o_bar_ref"] / f,
    }

def symbolic():
    import sympy as sp
    k_bar_d_ref, c_bar_0_ref, b_bar_ref, lambda_o_bar_ref, f_gamma = sp.symbols(
        "k_bar_d_ref c_bar_0_ref b_bar_ref lambda_o_bar_ref f_gamma", real=True
    )
    return {
        "k_bar_d": k_bar_d_ref / f_gamma,
        "c_bar_0": c_bar_0_ref / f_gamma,
        "b_bar": b_bar_ref / f_gamma,
        "lambda_o_bar": lambda_o_bar_ref / f_gamma,
    }
