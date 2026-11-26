import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _timeseries import read_timeseries_curve

def read(image, at):
    return read_timeseries_curve(image, at, swatch_idx=5, y_ticks_vals=[0.5, 0.4, 0.3, 0.2, 0.1, 0.0], rgb_tol=45.0)
