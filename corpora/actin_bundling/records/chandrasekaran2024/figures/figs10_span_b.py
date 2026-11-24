import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs5_10 import read_ellipsoid_boxplot
def read(image, at): return read_ellipsoid_boxplot(image, [1.5,1.0,0.5,0.0], ["Ring","Shell-SB"], "Supp. Fig. 10 span b")
