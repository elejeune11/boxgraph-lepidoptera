def run(inputs, settings):
    R = inputs["R_gamma"]
    return {
        "k_bar_d": inputs["k_bar_d_ref"] * R,
        "c_bar_0": inputs["c_bar_0_ref"] * R,
        "b_bar": inputs["b_bar_ref"] * R,
        "lambda_o_bar": inputs["lambda_o_bar_ref"] * R,
    }

def symbolic():
    import sympy as sp
    k_bar_d_ref, c_bar_0_ref, b_bar_ref, lambda_o_bar_ref, R_gamma = sp.symbols(
        "k_bar_d_ref c_bar_0_ref b_bar_ref lambda_o_bar_ref R_gamma", real=True
    )
    return {
        "k_bar_d": k_bar_d_ref * R_gamma,
        "c_bar_0": c_bar_0_ref * R_gamma,
        "b_bar": b_bar_ref * R_gamma,
        "lambda_o_bar": lambda_o_bar_ref * R_gamma,
    }
