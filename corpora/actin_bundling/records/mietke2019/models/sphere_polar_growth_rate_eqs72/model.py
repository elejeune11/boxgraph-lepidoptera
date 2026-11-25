def run(inputs, settings):
    return {
        "lambda_1": (
            -2.0 * inputs["D"] / (inputs["R0"] ** 2)
            - inputs["k"]
            + (inputs["xi"] / inputs["eta_b"]) * inputs["c0_df_dc0"]
        )
    }


def symbolic():
    import sympy as sp

    D, R0, k, xi, eta_b, c0_df_dc0 = sp.symbols("D R0 k xi eta_b c0_df_dc0", nonnegative=True)
    return {"lambda_1": -2 * D / R0**2 - k + (xi / eta_b) * c0_df_dc0}
