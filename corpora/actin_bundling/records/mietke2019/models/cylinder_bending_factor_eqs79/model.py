def run(inputs, settings):
    k_n = inputs["k_n"]
    r0 = inputs["r0"]
    kr2 = (k_n * r0) ** 2
    bc = 4.0 * (k_n ** 2) / (r0 ** 2) + ((2.0 * kr2 - 3.0) * (kr2 - 1.0)) / (r0 ** 4)
    return {"B_c": bc}


def symbolic():
    import sympy as sp

    k_n, r0 = sp.symbols("k_n r0", positive=True)
    return {"B_c": 4 * k_n**2 / r0**2 + (2 * k_n**2 * r0**2 - 3) * (k_n**2 * r0**2 - 1) / r0**4}
