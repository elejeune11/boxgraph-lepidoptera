import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import _fig1_2 as M
def read(image, at): return M.read_fig2d(image)
