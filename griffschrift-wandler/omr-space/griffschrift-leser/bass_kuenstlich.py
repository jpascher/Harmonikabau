# Künstliche Buchstaben A–H, a–h in verschiedenen Schriften als zusätzliche Lerndaten (Merkmale wie bassread.feat_box)
import sys, numpy as np, cv2, random
from PIL import Image, ImageDraw, ImageFont
sys.path.insert(0, '/home/claude/bass'); import bassread as BR
F = ['/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-regular.otf', '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyrepagella-regular.otf',
     '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreschola-regular.otf', '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyrebonum-regular.otf',
     '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-regular.otf', '/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
     '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf', '/usr/share/fonts/truetype/crosextra/Caladea-Regular.ttf', '/usr/share/fonts/truetype/freefont/FreeSans.ttf',
     '/usr/share/fonts/truetype/baskerville/GFSBaskerville.otf', '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyrecursor-regular.otf',
     '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bold.otf','/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyrebonum-bold.otf','/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreschola-bold.otf','/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyrepagella-bold.otf']
random.seed(1); d = 18.0; XS = []; YS = []
for fp in F:
    for L in 'ABCDEFGHabcdefgh':
        for k in range(20):
            cap = random.uniform(0.9, 2.4) * d; size = int(cap / 0.7)
            fnt = ImageFont.truetype(fp, size); im = Image.new('L', (size * 3, size * 3), 255); dr = ImageDraw.Draw(im)
            dr.text((size, size), L, font=fnt, fill=0)
            g = np.array(im)
            if random.random() < 0.4: sc=random.uniform(0.35,0.7); g = cv2.resize(cv2.resize(g,None,fx=sc,fy=sc,interpolation=cv2.INTER_AREA),(g.shape[1],g.shape[0]),interpolation=cv2.INTER_LINEAR)
            if random.random() < 0.5: g = cv2.GaussianBlur(g, (0, 0), random.uniform(0.3, 1.0))
            bw = (g < random.randint(120, 200)).astype(np.uint8) * 255
            if random.random() < 0.3: bw = cv2.dilate(bw, np.ones((2, 2), np.uint8))
            n, lab, S, _ = cv2.connectedComponentsWithStats(bw)
            if n < 2: continue
            i = 1 + int(np.argmax(S[1:, 4])); x, y, w, h, a = S[i]
            m = ((lab[y:y + h, x:x + w] == i) * 255).astype(np.uint8)
            ys4 = y - random.uniform(1.5, 7.5) * d
            XS.append(BR.feat_box(m, d, ys4, (x, y, w, h))); YS.append(L)
np.savez_compressed('/home/claude/bass/synth.npz', X=np.array(XS), Y=np.array(YS), K=np.array(['S'] * len(YS)))
print(len(YS))
