"""Extract a single SAXS peak wavevector by ordinal index from a digitized peak-list curve."""


def run(inputs, settings):
    curve = inputs["peaks"]
    target = int(settings["peak_index"])
    for x, y in zip(curve["x"], curve["y"]):
        if int(round(x)) == target:
            return {"q_peak": float(y)}
    raise ValueError(f"peak_index {target} not found in peaks x={curve['x']}")
