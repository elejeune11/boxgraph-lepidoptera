import importlib.util
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

_base = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("_sphere_arc_cphi", _base / "sphere_active_surface_arc_profiles" / "model.py")
_arc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_arc)


def run(inputs, settings):
    xs = [float(x) for x in inputs["sample_curve"]["x"]]
    res = _arc.solve_with_c_init(inputs, settings)
    ys = CubicSpline(res["s_norm"], res["Cphi"])(np.array(xs, dtype=float))
    return {"Cphi_profile": {"x": xs, "y": [float(v) for v in ys]}}
