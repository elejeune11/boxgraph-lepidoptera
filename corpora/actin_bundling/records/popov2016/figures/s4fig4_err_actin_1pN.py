import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _s4_excel import read_s4_excel

def read(image, at):
    return read_s4_excel(
        image, at, series='blue',
        y_grid_vals=[0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.0],
        x_grid_vals=[1.0, 10.0, 20.0, 40.0],
        legend_frac_box=(0.04, 0.10, 0.30, 0.32),
        has_error_bars=False,
        fallback_on_green_xs=(1.0,),
        non_negative=True,
    )
