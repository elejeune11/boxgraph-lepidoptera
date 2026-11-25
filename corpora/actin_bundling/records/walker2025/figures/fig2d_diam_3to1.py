import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _freq_poly import read_freq_curve

def read(image, at):
    return read_freq_curve(
        image,
        x_tick_vals=[0.5 * k for k in range(18)],
        y_tick_vals=[160.0, 120.0, 80.0, 40.0, 0.0],
        color_kind="blue",
        at=at,
        label="Fig. 2D 3:1 mini-Lpd:VASP",
    )
