def run(inputs, settings):
    k = inputs["k"]
    ell = inputs["ell"]
    denom = inputs["gamma"] * (1.0 + inputs["nu"] * (k * ell) ** 2)
    return {"v_hat_amp": k * inputs["zdm0"] * inputs["df_hat"] / denom}


def symbolic():
    import sympy as sp

    k, zdm0, gamma, nu, ell = sp.symbols("k zdm0 gamma nu ell", positive=True)
    df_hat = sp.symbols("df_hat", real=True)
    return {"v_hat_amp": k * zdm0 * df_hat / (gamma * (1 + nu * k**2 * ell**2))}
