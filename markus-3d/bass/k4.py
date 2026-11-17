import re
src=open('hebel.py').read()
cut=src.index('# Reihe 4: 11 Klappen')
exec_src=src[:cut]
