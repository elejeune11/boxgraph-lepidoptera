def run(inputs, settings):
    # w_bundle in um -> nm by multiplying by 1000; d_spacing in nm
    return {"N_fil": 1000.0 * inputs["w_bundle"] / inputs["d_spacing"]}


def symbolic():
    import sympy as sp

    w_bundle, d_spacing = sp.symbols("w_bundle d_spacing", positive=True)
    return {"N_fil": sp.Integer(1000) * w_bundle / d_spacing}
