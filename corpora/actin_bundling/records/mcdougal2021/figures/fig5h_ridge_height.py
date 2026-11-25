"""Calibrated reader for Fig. 5H top (p. 5): pooled mean +/- SD ridge peak height (nm) vs development (%)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5h import calibrate_fig5h, find_purple_stems, read_bars_from_stems


def read(image, at):
    x_axis, y_height, _y_spacing, purple, _teal, box, r_sep = calibrate_fig5h(image)
    stems = find_purple_stems(purple, box, r_sep)
    return read_bars_from_stems(purple, x_axis, y_height, stems, int(box[1]) + 60, int(r_sep) - 5)
