import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pca_sil import read_silhouette
def read(image, at): return read_silhouette(image, [0.25, 0.0], [2,4,6,8,10], "Supp. Fig. 9C")
