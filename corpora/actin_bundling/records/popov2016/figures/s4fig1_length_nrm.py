import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _s4_excel import read_s4_excel

def read(image, at):
    return read_s4_excel(
        image, at, series='blue',
        y_grid_vals=[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0],
        x_grid_vals=list(range(11)),
        legend_frac_box=(0.68, 0.35, 0.95, 0.72),
        has_error_bars=False,
        fallback_on_red_xs=(0.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0),
        non_negative=True,
    )
