import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import _fig5_S1 as M
def read(image, at): return M.read_tactoid_data(image, "figS1_right")
