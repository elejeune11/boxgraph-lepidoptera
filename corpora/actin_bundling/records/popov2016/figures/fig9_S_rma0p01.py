import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _heatmaps import read_heatmap_row

def read(image, at):
    return read_heatmap_row(image, at, row_idx=1, cb_labels=[0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1], min_half_width=0.008)
