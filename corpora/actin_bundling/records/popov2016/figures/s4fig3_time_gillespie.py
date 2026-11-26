import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _s4_excel import read_s4_excel

def read(image, at):
    return read_s4_excel(
        image, at, series='red',
        y_grid_vals=[7000.0, 6000.0, 5000.0, 4000.0, 3000.0, 2000.0, 1000.0, 0.0],
        x_grid_vals=list(range(11)),
        legend_frac_box=(0.68, 0.35, 0.95, 0.72),
        has_error_bars=True,
        non_negative=True,
    )
