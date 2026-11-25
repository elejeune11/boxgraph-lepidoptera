import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bar_replicates import read_bar_tops_and_replicates
def read(image, at):
    return read_bar_tops_and_replicates(image, y_tick_vals=[40.0, 35.0, 30.0, 25.0, 20.0, 15.0, 10.0, 5.0, 0.0], x_vals=[1.0, 2.0, 3.0], max_res_px=5.5, label='Fig. 7D VASP actin enrichment')
