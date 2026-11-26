def run(inputs, settings):
    return {
        "U_xl": (
            inputs["U_xl_stretch"]
            + inputs["U_xl_bind"]
            + inputs["U_xl_confine"]
            + inputs["U_xl_align"]
        )
    }


def symbolic():
    import sympy as sp

    U_xl_stretch, U_xl_bind, U_xl_confine, U_xl_align = sp.symbols(
        "U_xl_stretch U_xl_bind U_xl_confine U_xl_align", real=True
    )
    return {"U_xl": U_xl_stretch + U_xl_bind + U_xl_confine + U_xl_align}
