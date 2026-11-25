import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bar_replicates import read_bar_tops_and_replicates

def read(image, at):
    return read_bar_tops_and_replicates(
        image,
        y_tick_vals=[60.0, 50.0, 40.0, 30.0, 20.0, 10.0, 0.0],
        x_vals=[50.0, 100.0, 150.0, 500.0],
        label="Fig. 1H mini-Lpd partition",
    )
