def run(inputs, settings):
    return {"n_repeats": inputs["domain_in_wavelengths"] / inputs["wavelength_ratio"]}


def symbolic():
    import sympy as sp
    domain_in_wavelengths, wavelength_ratio = sp.symbols(
        "domain_in_wavelengths wavelength_ratio", positive=True
    )
    return {"n_repeats": domain_in_wavelengths / wavelength_ratio}
