import math


def run(inputs, settings):
    L_mem_0 = float(inputs["L_mem_0"])  # um
    int_c0 = float(inputs["int_c0"])    # s
    f_mem = float(inputs["f_mem"])      # pN
    delta = float(inputs["delta"])      # nm
    N = float(inputs["N"])              # dimensionless
    kBT = float(inputs["kBT"])          # pN*nm
    L_star = float(inputs["L_star"])    # um
    v_0 = float(inputs["v_0"])          # um/s

    l0 = L_mem_0 / L_star
    boltz = math.exp(-f_mem * delta / (N * kBT))
    l_t = math.sqrt(1.0 + l0 ** 2 + 2.0 * l0 + 2.0 * (v_0 / L_star) * boltz * int_c0) - 1.0
    return {"L_mem": L_star * l_t}


def symbolic():
    import sympy as sp

    L_mem_0, int_c0, f_mem, delta, N, kBT, L_star, v_0 = sp.symbols(
        "L_mem_0 int_c0 f_mem delta N kBT L_star v_0", positive=True
    )
    l0 = L_mem_0 / L_star
    l_t = sp.sqrt(1 + l0**2 + 2 * l0 + 2 * (v_0 / L_star) * sp.exp(-f_mem * delta / (N * kBT)) * int_c0) - 1
    return {"L_mem": L_star * l_t}
