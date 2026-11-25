import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _ar_boxplot import read_ar_boxplot
def read(image, at):
    return read_ar_boxplot(image, y_tick_vals=[float(k) for k in range(16, -1, -1)], x_vals=[0.0, 0.5, 1.0, 1.25], metric='median', n_cols=5, select_indices=[0,1,2,3], label='Fig. 6J Eps15-Lifeact AR')
