import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bar_replicates import read_bar_tops_and_replicates
def read(image, at):
    return read_bar_tops_and_replicates(image, y_tick_vals=[0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0], x_vals=[0.0, 1.0], max_res_px=5.5, label='Fig. 5J monomer-mini-Lpd frac AR>1.2')
