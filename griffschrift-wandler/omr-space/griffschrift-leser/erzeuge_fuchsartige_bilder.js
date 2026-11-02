// Übungsbilder im Fuchs-ähnlichen Stil: kleine dünne Kreuze, graue Linien, geringe Auflösung (Linienabstand ~6 px), zufällig variiert
const {chromium}=require('playwright');const fs=require('fs');
(async()=>{const b=await chromium.launch();
const demos=JSON.parse(fs.readFileSync('demos.json'));let k=0;
for(const dsf of [0.75,0.85,1.0]){
 const pg=await (await b.newContext({locale:'de-DE',viewport:{width:1100,height:1000},deviceScaleFactor:dsf,colorScheme:'light'})).newPage();
 await pg.goto('http://localhost:8765/griffschrift-wandler/index.html?lang=de');await pg.waitForTimeout(3000);
 for(const gs of ['ramoner','fuchs'])for(const d of demos){
  await pg.evaluate(async([d,gs])=>{const g=document.getElementById('gstyle');g.value=gs;g.dispatchEvent(new Event('change'));const fl=document.getElementById('fill');fl.value='0';fl.dispatchEvent(new Event('change'));await new Promise(r=>setTimeout(r,300));await loadDemo(d);},[d,gs]);
  await pg.waitForTimeout(1000);
  const n=await pg.evaluate(()=>document.querySelectorAll('#griff svg.gsys').length);
  for(let i=0;i<Math.min(n,4);i++){
   const gt=await pg.evaluate(([i,seed])=>{let s0=seed;const rnd=()=>{s0=(s0*9301+49297)%233280;return s0/233280;};
    const s=document.querySelectorAll('#griff svg.gsys')[i];
    // Stil verändern
    const cs=0.55+0.35*rnd(),sw=0.8+0.9*rnd(),gray=Math.round(60+120*rnd()),es=0.8+0.25*rnd();
    s.querySelectorAll('line.sl').forEach(l=>{l.style.stroke=`rgb(${gray},${gray},${gray})`;l.style.strokeWidth=0.6+0.6*rnd();});
    const heads=[];
    s.querySelectorAll('g.gn ellipse').forEach(e=>{const cx=+e.getAttribute('cx'),cy=+e.getAttribute('cy');e.setAttribute('rx',5.4*es);e.setAttribute('ry',3.9*es*(0.9+0.3*rnd()));heads.push({x:cx,y:cy,t:'o',fill:e.style.fill!=='none'});});
    s.querySelectorAll('g.gn path').forEach(p=>{const m=(p.getAttribute('d')||'').match(/^M([\d.]+) ([\d.]+) L([\d.]+) ([\d.]+) M/);if(!m)return;const cx=(+m[1]+ +m[3])/2,cy=(+m[2]+ +m[4])/2,r=4.5*cs;
      p.setAttribute('d',`M${cx-r} ${cy-r} L${cx+r} ${cy+r} M${cx-r} ${cy+r} L${cx+r} ${cy-r}`);p.setAttribute('stroke-width',sw);heads.push({x:cx,y:cy,t:'x'});});
    s.querySelectorAll('g.gn circle:not(.hh)').forEach(c=>{if(c.getAttribute('fill')==='none')c.setAttribute('r',5.6*cs);});
    s.querySelectorAll('text.mn,rect.brk,rect.msel,rect.bhit,rect.chit').forEach(t=>t.style.display=rnd()<0.5?'none':'');
    const bass=[...s.querySelectorAll('text.bl')].map(t=>({x:+t.getAttribute('x'),y:+t.getAttribute('y'),s:t.textContent}));
    const dl=[...s.querySelectorAll('line.dl')].map(l=>({x1:+l.getAttribute('x1'),x2:+l.getAttribute('x2'),y:+l.getAttribute('y1')}));
    const sl=[...s.querySelectorAll('line.sl')].filter(l=>l.getAttribute('y1')===l.getAttribute('y2')).map(l=>+l.getAttribute('y1'));
    const neg=[];s.querySelectorAll('text.mn,text.vl,text.ts,text.pw').forEach(t=>{if(t.style.display==='none')return;const bb=t.getBBox();neg.push({x:bb.x+bb.width/2,y:bb.y+bb.height/2});});
    s.querySelectorAll('circle.rp').forEach(c=>neg.push({x:+c.getAttribute('cx'),y:+c.getAttribute('cy')}));
    const byUid={};P.measures.forEach(m=>m.ev.forEach(e=>{if(e._uid)byUid[e._uid]=e;}));
    const grips=[...s.querySelectorAll('g.gn[data-uid]')].map(g=>{const e=byUid[g.dataset.uid];const r=e&&e._r;if(!r||!r.g)return null;const bb=g.getBBox();return {x:bb.x+bb.width/2,dir:r.g.dir,btns:r.g.btns.map(q=>[q.b.row,q.b.num]),pitches:r.g.btns.map(q=>q.m),dur:e.dur};}).filter(Boolean);
    return {heads,bass,dl,neg,grips,staff:[...new Set(sl)].sort((a,b)=>a-b)};},[i,k*7+13]);
   const el=(await pg.$$('#griff svg.gsys'))[i];const fn=`/tmp/claude-0/gs2/f${String(k).padStart(3,'0')}`;await el.screenshot({path:fn+'.png'});
   fs.writeFileSync(fn+'.json',JSON.stringify({dsf,...gt}));k++;}}
 await pg.close();}
console.log('Zeilen',k);await b.close();})();
