def run(inputs, settings):
    # k_bend in pN*um^2, kBT in pN*nm (= 1e-3 pN*um) -> L_p in um
    return {"L_p": float(inputs["k_bend"]) / (float(inputs["kBT"]) * 1e-3)}


def symbolic():
    import sympy as sp

    k_bend, kBT = sp.symbols("k_bend kBT", positive=True)
    return {"L_p": k_bend / (kBT * sp.Rational(1, 1000))}
