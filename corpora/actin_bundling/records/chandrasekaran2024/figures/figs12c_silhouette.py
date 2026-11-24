import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pca_sil import read_silhouette
def read(image, at): return read_silhouette(image, [1.0, 0.8, 0.6, 0.4], [2,3,4,5,6,7,8,9,10], "Supp. Fig. 12C")
