import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _freq_poly import read_freq_curve

def read(image, at):
    return read_freq_curve(
        image,
        x_tick_vals=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        y_tick_vals=[500.0, 400.0, 300.0, 200.0, 100.0, 0.0],
        color_kind="red",
        at=at,
        label="Fig. 1I 100 mM NaCl",
    )
