import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _ar_boxplot import read_ar_boxplot
def read(image, at):
    return read_ar_boxplot(image, y_tick_vals=[4.25 - 0.25*k for k in range(15)], x_vals=[2.0], metric='mean', n_cols=5, select_indices=[4], label='Fig. 4C mini-Lpd + LatA AR')
