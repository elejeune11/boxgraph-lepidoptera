import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _s4_excel import read_s4_excel

def read(image, at):
    return read_s4_excel(
        image, at, series='blue',
        y_grid_vals=[6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
        x_grid_vals=[1.0, 10.0, 20.0, 40.0],
        legend_frac_box=(0.55, 0.15, 0.95, 0.42),
        has_error_bars=False,
        non_negative=True,
    )
