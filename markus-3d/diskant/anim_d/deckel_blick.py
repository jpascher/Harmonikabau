import pickle, numpy as np, pyvista as pv
pv.global_theme.background="white"
D=pickle.load(open("netze.pkl","rb")); T=D["teile"]
def poly(v,t): return pv.PolyData(v, np.hstack([np.full((len(t),1),3),t]).astype(np.int64))
p=pv.Plotter(off_screen=True,window_size=(1800,1200),shape=(2,2))
cams=[[(190,-300,-60),(190,0,-60),(0,0,-1)],[(20,-120,-120),(55,0,-60),(0,0,-1)],[(360,-120,-120),(330,0,-60),(0,0,-1)],[(190,-150,-200),(190,-10,-100),(0,0,-1)]]
for i,cam in enumerate(cams):
    p.subplot(i//2,i%2)
    for n,(f,v,t) in T.items():
        if n.startswith("Griffstock_Koerper"): p.add_mesh(poly(v,t),color="burlywood")
        elif n.startswith("Griffstock_Abdeckung"): p.add_mesh(poly(v,t),color="steelblue",opacity=0.6)
        elif n.startswith(("Kunststoff","Dreieck")): p.add_mesh(poly(v,t),color="gray")
    p.camera_position=cam
p.screenshot("deckel_blick.png")
