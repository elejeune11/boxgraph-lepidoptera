import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bar_replicates import read_bar_tops_and_replicates
def read(image, at):
    return read_bar_tops_and_replicates(image, y_tick_vals=[1.0 - 0.1*k for k in range(11)], x_vals=[1.0, 2.0, 3.0], max_res_px=5.5, label='Fig. 7G VASP mutants frac AR>1.2')
