"""Panel crop definitions (at 1200 dpi) and printed tick-label specifications for tang1996.

Contains strictly the axis tick labels visible on each printed figure panel — no data values.
"""

from __future__ import annotations

import numpy as np
from _common import PanelCal, calibrate_panel

# (page, dpi, crop_at_1200dpi)
PANEL_SPECS: dict[str, dict] = {
    "fig1": dict(
        page=2, dpi=1200, crop=(5440, 940, 9080, 3600),
        x_scale="log10", y_scale="log10",
        x_major=[0.1, 1.0, 10.0],
        x_edges=(None, 100.0),
        x_all=[0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
               1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
               10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0],
        y_major=[1.0, 0.1],
        y_edges=(None, None),
        y_all=[2.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
               0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04],
    ),
    "fig2a": dict(
        page=3, dpi=1200, crop=(1120, 960, 4800, 3700),
        x_scale="linear", y_scale="log10",
        x_major=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        x_edges=(None, None),
        x_all=[round(-0.05 + 0.05 * k, 4) for k in range(23)],
        y_major=[1.0, 0.1],
        y_edges=(None, None),
        y_all=[3.0, 2.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
               0.1, 0.09, 0.08, 0.07, 0.06, 0.05],
    ),
    "fig2b": dict(
        page=3, dpi=1200, crop=(1120, 3900, 4800, 6700),
        x_scale="linear", y_scale="log10",
        x_major=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
        x_edges=(None, None),
        x_all=[float(2 * k) for k in range(26)],
        y_major=[1.0, 0.1],
        y_edges=(None, None),
        y_all=[4.0, 3.0, 2.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
               0.1, 0.09, 0.08, 0.07, 0.06],
    ),
    "fig3a": dict(
        page=3, dpi=1200, crop=(5300, 900, 9140, 3880),
        x_scale="linear", y_scale="log10",
        x_major=[0.0, 20.0, 40.0, 60.0, 80.0],
        x_edges=(None, 100.0),
        x_all=[float(-15 + 5 * k) for k in range(23)],
        x_sides=("top",),  # bottom axis has 3 dotted onset lines hitting it
        y_major=[1.0, 0.1],
        y_edges=(10.0, 0.01),
        y_all=[9, 8, 7, 6, 5, 4, 3, 2, 1.0,
               0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
               0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02],
    ),
    "fig3b": dict(
        page=3, dpi=1200, crop=(5300, 3900, 9140, 6860),
        x_scale="linear", y_scale="linear",
        x_major=[50.0, 100.0],
        x_edges=(0.0, 150.0),
        x_all=[float(5 * k) for k in range(1, 30)],
        y_major=[50.0, 40.0, 30.0, 20.0, 10.0],
        y_edges=(60.0, 0.0),
        y_all=[float(60 - 2 * k) for k in range(1, 30)],
    ),
    "fig4": dict(
        page=4, dpi=1200, crop=(5440, 920, 9140, 3660),
        x_scale="linear", y_scale="linear",
        x_major=[0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0],
        x_edges=(None, None),
        x_all=[float(-15 + 5 * k) for k in range(28)],
        x_sides=("top",),  # bottom tick near x=10 touches low data markers
        y_major=[1.0, 0.8, 0.6, 0.4, 0.2],
        y_edges=(1.2, 0.0),
        y_all=[round(1.1 - 0.1 * k, 2) for k in range(11)],
    ),
    "fig5a": dict(
        page=5, dpi=1200, crop=(1200, 920, 4680, 3620),
        x_scale="linear", y_scale="log10",
        x_major=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        x_edges=(None, None),
        x_all=[round(0.2 * k, 2) for k in range(30)],
        y_major=[0.1],
        y_edges=(1.0, None),
        y_all=[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06],
    ),
    "fig5b": dict(
        page=5, dpi=1200, crop=(1200, 3820, 4680, 6380),
        x_scale="log10", y_scale="log10",
        x_major=[100.0],
        x_edges=(None, 200.0),
        x_all=[30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
        x_match_by_all=True,
        y_major=[10.0, 1.0],
        y_edges=(None, None),
        y_all=[10.0, 9, 8, 7, 6, 5, 4, 3, 2, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5],
    ),
    "fig6a": dict(
        page=5, dpi=1200, crop=(5320, 920, 9100, 3740),
        x_scale="linear", y_scale="linear",
        x_major=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        x_edges=(None, None),
        x_all=[round(-0.3 + 0.1 * k, 2) for k in range(37)],
        y_major=[0.3, 0.2, 0.1],
        y_edges=(0.4, 0.0),
        y_all=[round(0.38 - 0.02 * k, 2) for k in range(19)],
    ),
    "fig6b": dict(
        page=5, dpi=1200, crop=(5320, 3900, 9100, 6740),
        x_scale="log10", y_scale="log10",
        x_major=[1.0, 10.0],
        x_edges=(0.1, 100.0),
        x_all=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
               2, 3, 4, 5, 6, 7, 8, 9, 10.0,
               20, 30, 40, 50, 60, 70, 80, 90],
        y_major=[0.1],
        y_edges=(1.0, None),
        y_all=[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06],
    ),
    "fig7": dict(
        page=5, dpi=1200, crop=(5440, 7720, 9080, 10400),
        x_scale="log10", y_scale="log10",
        x_major=[0.001, 0.01],
        x_edges=(0.0001, 0.1),
        x_all=[2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3,
               2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2,
               2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2],
        y_major=[1.0, 0.1],
        y_edges=(None, 0.01),
        y_all=[4.0, 3.0, 2.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
               0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02],
        y_sides=("right",),  # left curve starts right on left border near y=0.11
    ),
    "fig8a": dict(
        page=6, dpi=1200, crop=(1120, 940, 4800, 3720),
        x_scale="linear", y_scale="log10",
        x_major=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        x_edges=(None, 6.0),
        x_all=[round(0.2 * k, 2) for k in range(30)],
        y_major=[1.0, 0.1],
        y_edges=(None, None),
        y_all=[1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
               0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04],
    ),
    "fig8b": dict(
        page=6, dpi=1200, crop=(1120, 3960, 4800, 6720),
        x_scale="log10", y_scale="log10",
        x_major=[1.0, 10.0],
        x_edges=(None, None),
        x_all=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10.0, 20.0],
        y_major=[1.0, 0.1],
        y_edges=(None, None),
        y_all=[1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
               0.1, 0.09, 0.08, 0.07, 0.06, 0.05],
    ),
}


def calibrate(name: str, image: np.ndarray) -> PanelCal:
    s = PANEL_SPECS[name]
    kw = {k: v for k, v in s.items() if k not in ("page", "dpi", "crop")}
    return calibrate_panel(image, **kw)
