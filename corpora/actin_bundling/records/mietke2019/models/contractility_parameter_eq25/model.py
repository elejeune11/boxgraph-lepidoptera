def run(inputs, settings):
    return {"alpha": inputs["xi_over_gamma"] * inputs["c0_df_dc0"]}


def symbolic():
    import sympy as sp

    xi_over_gamma, c0_df_dc0 = sp.symbols("xi_over_gamma c0_df_dc0", real=True)
    return {"alpha": xi_over_gamma * c0_df_dc0}
