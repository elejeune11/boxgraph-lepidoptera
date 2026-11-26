import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _timeseries import read_timeseries_curve

def read(image, at):
    return read_timeseries_curve(image, at, swatch_idx=5, y_ticks_vals=[1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6], rgb_tol=45.0)
