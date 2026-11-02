// Griffschrift-Zeilen als PNG + Lösung (Köpfe, Bässe, Druckstriche, Taktstriche) aus dem SVG
const {chromium}=require('playwright');const fs=require('fs');
(async()=>{const b=await chromium.launch();const pg=await (await b.newContext({locale:'de-DE',viewport:{width:1400,height:1000},deviceScaleFactor:2,colorScheme:'light'})).newPage();
await pg.goto('http://localhost:8765/griffschrift-wandler/index.html?lang=de');await pg.waitForTimeout(3000);
const demos=JSON.parse(fs.readFileSync('demos.json'));let k=0;
for(const gs of ['ramoner','fuchs'])for(const d of demos){await pg.evaluate(async([d,gs])=>{const g=document.getElementById('gstyle');g.value=gs;g.dispatchEvent(new Event('change'));const fl=document.getElementById('fill');fl.value='0';fl.dispatchEvent(new Event('change'));await new Promise(r=>setTimeout(r,300));await loadDemo(d);},[d,gs]);await pg.waitForTimeout(1200);
 const n=await pg.evaluate(()=>document.querySelectorAll('#griff svg.gsys').length);
 for(let i=0;i<Math.min(n,6);i++){
  const gt=await pg.evaluate(i=>{const s=document.querySelectorAll('#griff svg.gsys')[i];const heads=[];
   s.querySelectorAll('g.gn ellipse').forEach(e=>heads.push({x:+e.getAttribute('cx'),y:+e.getAttribute('cy'),t:'o',fill:e.style.fill!=='none'}));
   s.querySelectorAll('g.gn path').forEach(p=>{const m=(p.getAttribute('d')||'').match(/^M([\d.]+) ([\d.]+) L([\d.]+) ([\d.]+) M/);if(m)heads.push({x:(+m[1]+ +m[3])/2,y:(+m[2]+ +m[4])/2,t:'x'});});
   const bass=[...s.querySelectorAll('text.bl')].map(t=>({x:+t.getAttribute('x'),y:+t.getAttribute('y'),s:t.textContent}));
   const dl=[...s.querySelectorAll('line.dl')].map(l=>({x1:+l.getAttribute('x1'),x2:+l.getAttribute('x2'),y:+l.getAttribute('y1')}));
   const sl=[...s.querySelectorAll('line.sl')].filter(l=>l.getAttribute('y1')===l.getAttribute('y2')).map(l=>+l.getAttribute('y1'));
   const neg=[];s.querySelectorAll('text.mn,text.vl,text.ts,text.pw').forEach(t=>{const bb=t.getBBox();neg.push({x:bb.x+bb.width/2,y:bb.y+bb.height/2});});
   s.querySelectorAll('circle.rp').forEach(c=>neg.push({x:+c.getAttribute('cx'),y:+c.getAttribute('cy')}));
   s.querySelectorAll('g.gn').forEach(g=>{if(!g.querySelector('ellipse')&&!g.querySelector('path[d^="M"]')){const bb=g.getBBox();if(bb.width)neg.push({x:bb.x+bb.width/2,y:bb.y+bb.height/2});}});
   const byUid={};P.measures.forEach(m=>m.ev.forEach(e=>{if(e._uid)byUid[e._uid]=e;}));
   const grips=[...s.querySelectorAll('g.gn[data-uid]')].map(g=>{const e=byUid[g.dataset.uid];const r=e&&e._r;if(!r||!r.g)return null;const bb=g.getBBox();
     return {x:bb.x+bb.width/2,dir:r.g.dir,btns:r.g.btns.map(q=>[q.b.row,q.b.num]),pitches:r.g.btns.map(q=>q.m),dur:e.dur};}).filter(Boolean);
   return {heads,bass,dl,neg,grips,staff:[...new Set(sl)].sort((a,b)=>a-b)};},i);
  const el=(await pg.$$('#griff svg.gsys'))[i];const fn=`/tmp/claude-0/gs/z${String(k).padStart(3,'0')}`;await el.screenshot({path:fn+'.png'});
  fs.writeFileSync(fn+'.json',JSON.stringify({demo:d,...gt}));k++;}}
console.log('Zeilen',k);await b.close();})();
