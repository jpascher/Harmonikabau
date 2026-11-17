import math
exec(open('hebel_reg.py').read().split('\nif __name__')[0])
src=open('hebel.py').read()
i=src.index('_tmpdoc = ezdxf.new()'); j=src.index('\n# --- Klappenmitten')
exec(src[i:j])
L=math.hypot(*KNOPF_RICHT); d=(KNOPF_RICHT[0]/L,KNOPF_RICHT[1]/L); e=(d[1],-d[0])
ab=lambda p:(p[0]*d[0]+p[1]*d[1],p[0]*e[0]+p[1]*e[1])
kn=sorted(KNOPF,key=lambda q:-q[0])
A0=ab(kn[0])[0]-3.5; A={r:A0+4*(r-1) for r in (1,2,3,4)}; B={r:ab(kn[r-1])[1] for r in (1,2,3,4)}
print('A',A,'B',B)
for r in (1,2,3,4):
    w=[ab(p) for p in trafo(r,holz_ring(r))]
    print('Reihe',r,'Holz a max',round(max(q[0] for q in w),2),'bei b',round(max(w)[1],1), ' b-Bereich',round(min(q[1] for q in w),1),round(max(q[1] for q in w),1))
    for rr in (4,3,2,1):
        lo=B[rr+1]/2+B[rr]/2 if rr<4 else B[4]-30; hi=(B[rr]+B[rr-1])/2 if rr>1 else B[1]+9
        qs=[q[0] for q in w if lo<=q[1]<=hi]
        if qs: print('   ueber Stufe',rr,'(a=%.1f) Holz bis a=%.2f'%(A[rr],max(qs)))
