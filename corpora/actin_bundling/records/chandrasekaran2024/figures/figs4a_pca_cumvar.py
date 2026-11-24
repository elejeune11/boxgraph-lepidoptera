import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pca_sil import read_pca_cumvar
def read(image, at): return read_pca_cumvar(image, 5, "Supp. Fig. 4A")
