import pyvista as pv, sys
col = {"Holz": "burlywood", "Alu": "silver", "Verdeck": "#333333", "Leder": "saddlebrown", "Griffstock": "sienna", "Stahl": "steelblue", "Schraube": "goldenrod", "Klappe": "#4d3320", "Hebelholz": "peru", "Kunststoff": "ivory", "Knopf": "white"}
def bild(name, cam, clip=None):
    p = pv.Plotter(off_screen=True, window_size=(1600, 1000))
    for f, c in col.items():
        m = pv.read(f"vor_{f}.stl")
        if clip: m = m.clip(normal=(1,0,0), origin=(clip,0,0))
        p.add_mesh(m, color=c)
    p.set_background("white")
    if cam == "yz": p.view_yz()
    else: p.camera_position = cam
    p.screenshot(name)
bild(sys.argv[1]+"_a.png", [(560, -330, -230), (190, 100, 0), (0, 0, 1)])
bild(sys.argv[1]+"_b.png", [(190, -330, 330), (190, 100, 20), (0, 0, 1)])
bild(sys.argv[1]+"_s.png", "yz", clip=102)
