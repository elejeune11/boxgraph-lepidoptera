import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _logx_errbars import read_logx_red_markers

def read(image, at):
    return read_logx_red_markers(image, at, y_ticks_vals=[0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
