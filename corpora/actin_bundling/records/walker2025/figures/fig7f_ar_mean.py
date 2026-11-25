import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _ar_boxplot import read_ar_boxplot
def read(image, at):
    return read_ar_boxplot(image, y_tick_vals=[float(k) for k in range(16, -1, -1)], x_vals=[1.0, 2.0, 3.0], metric='mean', label='Fig. 7F VASP mutants AR')
