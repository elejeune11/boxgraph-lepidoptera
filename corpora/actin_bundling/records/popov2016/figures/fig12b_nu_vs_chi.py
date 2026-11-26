import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _logx_errbars import read_logx_red_markers

def read(image, at):
    return read_logx_red_markers(image, at, y_ticks_vals=[1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9])
