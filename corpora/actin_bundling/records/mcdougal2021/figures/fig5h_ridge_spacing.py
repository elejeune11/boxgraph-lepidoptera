"""Calibrated reader for Fig. 5H bottom (p. 5): mean +/- SD ridge periodicity/spacing (um) vs development (%)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5h import calibrate_fig5h, find_purple_stems, read_bars_from_stems


def read(image, at):
    x_axis, _y_height, y_spacing, purple, teal, box, r_sep = calibrate_fig5h(image)
    stems = find_purple_stems(purple, box, r_sep)
    return read_bars_from_stems(teal, x_axis, y_spacing, stems, int(r_sep) + 25, int(box[3]) - 10)
