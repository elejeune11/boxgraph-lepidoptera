def run(inputs, settings):
    xs = inputs["spacings"]["x"]
    ys = inputs["spacings"]["y"]
    by_idx = {int(round(x)): float(y) for x, y in zip(xs, ys)}
    return {
        "L_2A": by_idx[1],
        "L_2B": by_idx[2],
        "L_1A": by_idx[3],
        "L_1B": by_idx[4],
        "L_Vcardui": by_idx[5],
    }
