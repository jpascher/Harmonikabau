import sys, pickle, numpy as np
sys.argv = ["kammern.py", "boden2.dxf", "/tmp/claude-0/-home-claude/cf042130-e5d9-5fc2-925a-aa4079250231/scratchpad/anim/tmp2.step"]
G = {"__name__": "__main__"}
exec(open("kammern.py", encoding="utf-8").read(), G)
def mesh(s, tol=0.25):
    v, t = s.tessellate(tol, 0.4)
    return np.array([(p.x, p.y, p.z) for p in v], dtype=float), np.array(t, dtype=np.int64)
teile = [("Front", G["front"]), ("Rueckteil", G["rueck"]), ("Seitenteil_links", G["sl"]), ("Seitenteil_rechts", G["sr"]), ("Bassboden", G["boden"])]
teile += G["teile"] + G["mech"] + G["klappen"] + G["hoerner"] + G["fuesse"] + [("Frontblende", G["bl"].val())] + G["befest"]
alle = {}
for n, s in teile:
    try: alle[n] = mesh(s, 0.12 if s.BoundingBox().DiagonalLength < 40 else 0.3)
    except Exception as e: print("skip", n, e)
pickle.dump((alle, G["L"]), open("anim/explo.pkl", "wb")); print("Teile", len(alle))
