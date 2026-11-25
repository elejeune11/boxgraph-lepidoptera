import importlib.util
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

_base = Path(__file__).resolve().parent.parent
_arc_spec = importlib.util.spec_from_file_location("_sphere_arc_refp", _base / "sphere_active_surface_arc_profiles" / "model.py")
_arc = importlib.util.module_from_spec(_arc_spec)
_arc_spec.loader.exec_module(_arc)


def run(inputs, settings):
    xs = [float(x) for x in inputs["sample_curve"]["x"]]
    res = _arc.solve_with_c_init(inputs, settings)
    u_norm = res["u_norm"]
    x_arr = np.array(xs, dtype=float)
    return {
        "h_profile": {"x": xs, "y": [float(v) for v in CubicSpline(u_norm, res["h"])(x_arr)]},
        "Css_profile": {"x": xs, "y": [float(v) for v in CubicSpline(u_norm, res["Css"])(x_arr)]},
    }
