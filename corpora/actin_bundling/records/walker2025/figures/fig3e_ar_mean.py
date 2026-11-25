import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _ar_boxplot import read_ar_boxplot

def read(image, at):
    return read_ar_boxplot(
        image,
        y_tick_vals=[float(k) for k in range(11, -1, -1)],
        x_vals=[0.0, 0.25, 0.5, 0.75],
        metric="mean",
        label="Fig. 3E mini-Lpd+VASP AR",
    )
