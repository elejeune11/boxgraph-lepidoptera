import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bar_replicates import read_bar_tops_and_replicates

def read(image, at):
    return read_bar_tops_and_replicates(
        image,
        y_tick_vals=[1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0],
        x_vals=[0.0, 0.25, 0.5, 0.75],
        label="Fig. 3F fraction AR>1.2",
    )
