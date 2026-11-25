# 23.11.2026: Stub - der eigentliche Code steht in Pruefung\pruefen_code.py und wird bei jedem Lauf frisch gelesen
def run(context):
    g = {"__name__": "pruefen_code"}
    exec(open(r"C:\ACAD-DWG\0-Markus-3D\Pruefung\pruefen_code.py", encoding="utf-8").read(), g)
    g["run"](context)
