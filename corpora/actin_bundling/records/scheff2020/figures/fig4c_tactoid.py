import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import _fig3_4 as M
def read(image, at): return M.read_fig4c(image, "tactoid")
