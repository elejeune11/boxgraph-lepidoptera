import math


def run(inputs, settings):
    L_mem_0 = float(inputs["L_mem_0"])  # dimensionless as printed
    int_c0 = float(inputs["int_c0"])    # dimensionless as printed
    f_mem = float(inputs["f_mem"])      # pN
    delta = float(inputs["delta"])      # nm
    N = float(inputs["N"])              # dimensionless
    kBT = float(inputs["kBT"])          # pN*nm

    boltz = math.exp(-f_mem * delta / (N * kBT))
    L_mem = math.sqrt(1.0 + L_mem_0 ** 2 + 2.0 * L_mem_0 + 2.0 * boltz * int_c0) - 1.0
    return {"L_mem": L_mem}


def symbolic():
    import sympy as sp

    L_mem_0, int_c0, f_mem, delta, N, kBT = sp.symbols("L_mem_0 int_c0 f_mem delta N kBT", positive=True)
    return {
        "L_mem": sp.sqrt(1 + L_mem_0**2 + 2 * L_mem_0 + 2 * sp.exp(-f_mem * delta / (N * kBT)) * int_c0) - 1
    }
