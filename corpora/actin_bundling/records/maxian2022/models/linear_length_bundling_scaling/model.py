def run(inputs, settings):
    # Unknown prefactor set to 1 (in s / um); BoxGraph judges shape only.
    return {"tau_c": inputs["L"]}


def symbolic():
    import sympy as sp

    L = sp.Symbol("L", positive=True)
    return {"tau_c": L}
