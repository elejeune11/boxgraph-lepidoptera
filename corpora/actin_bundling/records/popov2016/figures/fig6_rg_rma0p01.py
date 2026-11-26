import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _heatmaps import read_heatmap_row

def read(image, at):
    return read_heatmap_row(image, at, row_idx=1, cb_labels=[1.8, 1.6, 1.4, 1.2, 1.0], min_half_width=0.012)
