import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _timeseries import read_timeseries_curve

def read(image, at):
    return read_timeseries_curve(image, at, swatch_idx=3, y_ticks_vals=[1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75], rgb_tol=45.0)
