def run(inputs, settings):
    ar_c = inputs["aspect_ratio_curve"]
    om_c = inputs["omega_target_curve"]
    xs_out = []
    ys_out = []
    for x_tgt in om_c["x"]:
        # Pair each plotted omega point with the nearest measured L/r point at the same filamin %
        best_i = min(range(len(ar_c["x"])), key=lambda i: abs(ar_c["x"][i] - x_tgt))
        if abs(ar_c["x"][best_i] - x_tgt) > 0.25:
            raise ValueError(f"No matching aspect-ratio measurement near filamin={x_tgt}%")
        ar = float(ar_c["y"][best_i])
        xs_out.append(float(x_tgt))
        ys_out.append(0.25 * ar * ar)
    return {"omega_curve": {"x": xs_out, "y": ys_out}}
