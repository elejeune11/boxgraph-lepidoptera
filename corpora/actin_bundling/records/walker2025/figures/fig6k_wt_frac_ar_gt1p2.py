import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bar_replicates import read_bar_tops_and_replicates
def read(image, at):
    return read_bar_tops_and_replicates(image, y_tick_vals=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0], x_vals=[1.0], n_bars=5, select_indices=[4], y_ticks_bottom_n=6, max_res_px=7.5, y_max_search=1.1, label='Fig. 6K WT Eps15 control frac AR>1.2')
