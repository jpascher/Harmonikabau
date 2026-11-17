import pickle, numpy as np, pyvista as pv
pv.global_theme.background="white"
D=pickle.load(open("netze.pkl","rb")); T=D["teile"]
def poly(v,t): return pv.PolyData(v, np.hstack([np.full((len(t),1),3),t]).astype(np.int64))
p=pv.Plotter(off_screen=True,window_size=(1800,1100),shape=(2,1))
for i,cam in enumerate([[(60,180,-170),(120,20,-70),(0,0,-1)],[(330,140,-60),(150,15,-70),(0,0,-1)]]):
    p.subplot(i,0)
    for n,(f,v,t) in T.items():
        if n.startswith("Griffstock_Koerper"): p.add_mesh(poly(v,t),color="burlywood")
        elif "/Holzteil" in n: p.add_mesh(poly(v,t),color="green" if "_R3_" in n else "orange" if "_R4_" in n else "red" if "_R1_" in n else "blue")
    p.camera_position=cam
p.screenshot("gs_blick.png")
