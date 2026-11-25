# 23.11.2026: Bassteil auf groessere Laenge strecken (BASS_DL = Mehrlaenge in mm, je Haelfte an beiden Enden).
# Wird von kammern.py vor dem Zusammenbau eingebunden; jedes Teil, das a.add() bekommt, wird gestreckt/verschoben:
# Ebenen X = STRECK_RAND bzw. L - STRECK_RAND (knapp innerhalb der Waende, dort laufen nur Bretter/Leisten durch);
# kurze Teile (< 30 mm) werden als Ganzes verschoben, lange an der Ebene geteilt und mit ihrem Querschnitt gefuellt.
import os as _os
from OCP.gp import gp_GTrsf, gp_Mat, gp_XYZ
from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform
BASS_DL = float(_os.environ.get("BASS_DL", "0") or 0)
STRECK_RAND = float(_os.environ.get("BASS_STRECK_RAND", "8.0"))
_streck_warn = []
def _halb(s_, x0_, x1_):
    bb_ = s_.BoundingBox()
    return s_.intersect(cq.Solid.makeBox(x1_ - x0_, bb_.ylen + 20, bb_.zlen + 20, cq.Vector(x0_, bb_.ymin - 10, bb_.zmin - 10)))
def _streck_ebene(s_, xc_, d_, name_):
    bb_ = s_.BoundingBox()
    if bb_.xmax <= xc_ + 1e-6: return s_
    if bb_.xmin >= xc_ - 1e-6: return s_.translate(cq.Vector(d_, 0, 0))
    if bb_.xlen < 30.0 or name_.startswith("Horn"):       # kleines Teil / Horn: als Ganzes nach seiner Mitte
        return s_.translate(cq.Vector(d_, 0, 0)) if (bb_.xmin + bb_.xmax) / 2 > xc_ else s_
    a_ = _halb(s_, xc_ - 1.0, xc_); b_ = _halb(s_, xc_, xc_ + 1.0)
    va_, vb_ = a_.Volume(), b_.Volume()
    if abs(va_ - vb_) > 0.02 * max(va_, vb_, 1e-6):
        _streck_warn.append((name_, round(xc_, 1), round(va_, 2), round(vb_, 2)))
    li_ = _halb(s_, bb_.xmin - 1, xc_); re_ = _halb(s_, xc_, bb_.xmax + 1).translate(cq.Vector(d_, 0, 0))
    g_ = gp_GTrsf(); g_.SetVectorialPart(gp_Mat(d_, 0, 0, 0, 1, 0, 0, 0, 1)); g_.SetTranslationPart(gp_XYZ(xc_ - (xc_ - 1.0) * d_, 0, 0))
    mi_ = cq.Shape.cast(BRepBuilderAPI_GTransform(a_.wrapped, g_, True).Shape())
    try: return li_.fuse(mi_, re_).clean()
    except Exception: return li_.fuse(mi_).fuse(re_)
# Front/Rueckteil tragen an den Enden Zinken, Ornamente (bis 29 mm) und Hornbohrungen (bis 45 mm): dort zwischen
# Hornbohrung und Knopfbereich strecken (glatter Bereich X 46..85); die Hoerner selbst wandern mit der Front.
FRONT_EBENE = float(_os.environ.get("BASS_FRONT_EBENE", "50.0"))   # zwischen Hornbohrung 1 (bis 45,3) und 2 (ab 55,8)
def strecke(obj_, name_=""):
    if BASS_DL <= 0 or not isinstance(obj_, cq.Shape): return obj_
    d2_ = BASS_DL / 2
    rand_ = FRONT_EBENE if name_.startswith(("Front", "Rueckteil", "Horn")) and not name_.startswith(("Frontleiste", "Front_Aufdoppelung", "Frontblende")) else STRECK_RAND
    s_ = _streck_ebene(obj_, L - rand_, d2_, name_)
    return _streck_ebene(s_, rand_, d2_, name_)
