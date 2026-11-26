"""Extract a single SAXS peak wavevector by ordinal peak_index from a digitized peak list."""


def run(inputs, settings):
    idx = int(settings.get("peak_index", 1))
    xs = inputs["peaks"]["x"]
    ys = inputs["peaks"]["y"]
    by_idx = {int(round(x)): float(y) for x, y in zip(xs, ys)}
    if idx not in by_idx:
        raise ValueError(f"peak_index {idx} not present in digitized peaks (available: {sorted(by_idx)})")
    return {"q_peak": by_idx[idx]}
