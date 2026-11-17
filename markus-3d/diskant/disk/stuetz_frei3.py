import numpy as np, pickle, collections
from shapely.geometry import box, Point
O=pickle.load(open("stuetz_obst2.pkl","rb")); D=pickle.load(open("stuetz_deck.pkl","rb"))
by=collections.defaultdict(list)
for x in np.arange(12,368,0.5):
    for y0 in np.arange(36,150,0.5):
        if not O.intersects(box(x-2.5,y0,x+2.5,y0+40)):
            # Abstand Mitte zum naechsten Lochrand
            m=Point(x,y0+20); inloch=not D.contains(m)
            by[x].append((y0,inloch, round(D.boundary.distance(m),1)))
for x in sorted(by):
    ys=by[x]; l=[y for y in ys if y[1]]
    print(x, ys[0][0], ys[-1][0], len(ys), "Mitte im Loch:", [(y[0],y[2]) for y in l][:6])
