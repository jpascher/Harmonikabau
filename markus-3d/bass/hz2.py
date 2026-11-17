import math
src=open('hebel.py').read()
# Kontext aus diskant.py-Konstanten
H_VD, T_VD = 32.0, 4.0
exec(src[:src.index('# --- Klappenmitten')])
i=src.index('knopf_u = sorted(KNOPF'); j=src.index('DREHPUNKTE = {}')
exec(src[i:j])
from shapely.geometry import Polygon as _SPG, LineString as _LS
for r in (1,2,3,4):
    g=holz_geclippt(r, holz_ring(r))
    out=[]
    for b in range(20,110,6):
        ln=_LS([_uy(-50,b),_uy(80,b)]); s=g.intersection(ln)
        if s.is_empty: continue
        amax=max(_ab(c)[0] for c in (s.coords if hasattr(s,'coords') else [c for q in s.geoms for c in q.coords]))
        fl=FRONT_CLIP.exterior.intersection(ln)
        af=max(_ab(c)[0] for c in (fl.coords if hasattr(fl,'coords') else [c for q in fl.geoms for c in q.coords])) if not fl.is_empty else None
        out.append((b, round(amax,1), round(af,1) if af else None))
    print(r, out)
