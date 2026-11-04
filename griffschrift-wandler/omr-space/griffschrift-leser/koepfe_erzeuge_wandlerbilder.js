// Griffschrift im Stil von Elke Margetich (4. Reihe als Dreieck, Spitze oben, hohl bei halben/ganzen Noten) als PNG + Lösung
// Aufruf: node gen_tri.js OUTDIR SEED ANTEIL piece1 piece2 ...   (piece = phx:… oder vm:… oder mid:…)
const {chromium}=require('playwright');const fs=require('fs');
const [out,seed0,N]=process.argv.slice(2);fs.mkdirSync(out,{recursive:true});const VM=fs.readdirSync('/tmp/claude-0/srv/midi/vm').filter(f=>f.endsWith('.mid')).map(f=>'vm:'+f.slice(0,-4));
(async()=>{const b=await chromium.launch();let k=0;let seed=+seed0;
 const rnd=()=>{seed=(seed*9301+49297)%233280;return seed/233280;};
 for(const dsf of [1,1.5,2]){
 const pg=await (await b.newContext({locale:'de-DE',viewport:{width:1100+Math.round(rnd()*300),height:1400},deviceScaleFactor:dsf,colorScheme:'light'})).newPage();
 pg.on('dialog',d=>d.dismiss());
 await pg.goto('http://127.0.0.1:8798/index.html?lang=de');await pg.waitForTimeout(2500);
 for(let ii=0;ii<+N;ii++){const pc=VM[Math.floor(rnd()*VM.length)];
  const vo=['2','3','4'][Math.floor(rnd()*3)];
  await pg.evaluate(async([pc,vo])=>{const set=(id,v)=>{const s=document.getElementById(id);s.value=v;s.dispatchEvent(new Event('change',{bubbles:true}));};
    {const st=document.getElementById('gstyle');const o=[...st.options].map(x=>x.value);set('gstyle',o[Math.floor(Math.random()*o.length)]);}set('fill','0');set('voices',vo);set('view','griff');
    const [kind,key]=pc.split(':');if(kind==='vm')await loadVM(key);else{const ds=document.getElementById('demo');ds.value=pc;await loadDemo(pc);}},[pc,vo]);
  await pg.waitForTimeout(1800);
  const n=await pg.evaluate(()=>document.querySelectorAll('#griff svg.gsys').length);
  for(let i=0;i<n;i++){
   const gt=await pg.evaluate(([i,sd])=>{let s0=sd;const rnd=()=>{s0=(s0*9301+49297)%233280;return s0/233280;};
    const s=document.querySelectorAll('#griff svg.gsys')[i];
    const sl=[...s.querySelectorAll('line.sl')].filter(l=>l.getAttribute('y1')===l.getAttribute('y2')).map(l=>+l.getAttribute('y1'));
    const ys=[...new Set(sl)].sort((a,b)=>a-b);if(ys.length<5)return null;const mid=ys[2];
    const byUid={};P.measures.forEach(m=>m.ev.forEach(e=>{if(e._uid)byUid[e._uid]=e;}));
    const gray=Math.round(40+150*rnd());s.querySelectorAll('line.sl').forEach(l=>{l.style.stroke=`rgb(${gray},${gray},${gray})`;l.style.strokeWidth=0.6+0.8*rnd();});
    const tw=9+4*rnd(),th=7.5+2.5*rnd(),tsw=0.9+0.9*rnd(),es=0.85+0.25*rnd(),xs=0.8+0.3*rnd(),xsw=1.2+1.2*rnd();
    s.querySelectorAll('text.mn,rect.brk,rect.msel,rect.bhit,rect.chit').forEach(t=>t.style.display=rnd()<0.5?'none':'');
    document.querySelectorAll('body *').forEach(el=>{if(getComputedStyle(el).position==='fixed')el.style.display='none';});
    s.querySelectorAll('[style]').forEach(el=>{const st=el.getAttribute('style');if(/var\(--(extra|muted|accent|zug|druck)\)/.test(st))el.setAttribute('style',st.replace(/var\(--(extra|muted|accent|zug|druck)\)/g,'#000'));});
    s.querySelectorAll('.cur,.selx,.forced').forEach(el=>el.classList.remove('cur','selx','forced'));
    const heads=[];let ntri=0;
    s.querySelectorAll('g.gn').forEach(g=>{const e=byUid[g.dataset.uid]||null;const src=e&&e.type==='tie'?e.src:e;const r=src&&src._r;
      const row4=new Set();if(r&&r.g)r.g.btns.forEach(q=>{if(q.b.row===4)row4.add(Math.round(mid-posOf(q.b)*4));});
      [...g.querySelectorAll('ellipse')].forEach(el=>{const cx=+el.getAttribute('cx'),cy=+el.getAttribute('cy');el.setAttribute('rx',5.4*es);el.setAttribute('ry',3.9*es);const st=el.getAttribute('style')||'';heads.push({x:cx,y:cy,t:'o',fill:!/fill:none/.test(st)});});
      [...g.querySelectorAll('path')].forEach(p=>{const m=(p.getAttribute('d')||'').match(/^M([-\d.]+) ([-\d.]+) L([-\d.]+) ([-\d.]+) M/);if(!m)return;
        const cx=(+m[1]+ +m[3])/2,cy=(+m[2]+ +m[4])/2,hollow=+(p.getAttribute('stroke-width')||2)<2;
        if(false){
          const circ=[...g.querySelectorAll('circle')].find(c=>!c.classList.contains('hh')&&Math.abs(+c.getAttribute('cx')-cx)<0.5&&Math.abs(+c.getAttribute('cy')-cy)<0.5);if(circ)circ.parentNode.removeChild(circ);
          const poly=document.createElementNS('http://www.w3.org/2000/svg','polygon');
          poly.setAttribute('points',`${cx-tw/2},${cy+th/2} ${cx+tw/2},${cy+th/2} ${cx},${cy-th/2}`);
          poly.setAttribute('style',hollow?`fill:none;stroke:#000;stroke-width:${tsw}`:'fill:#000;stroke:#000;stroke-width:0.5');
          p.replaceWith(poly);heads.push({x:cx,y:cy,t:'t',fill:!hollow});ntri++;}
        else{const rr=4.5*xs;p.setAttribute('d',`M${cx-rr} ${cy-rr} L${cx+rr} ${cy+rr} M${cx-rr} ${cy+rr} L${cx+rr} ${cy-rr}`);p.setAttribute('stroke-width',hollow?1.3:xsw);heads.push({x:cx,y:cy,t:'x',fill:!hollow});}});});
    return {heads,ntri,staff:ys};},[i,Math.floor(rnd()*100000)]);
   if(!gt||!gt.heads.length)continue;
   const el=(await pg.$$('#griff svg.gsys'))[i];const fn=`${out}/t${seed0}_${String(k).padStart(4,'0')}`;await el.screenshot({path:fn+'.png'});
   fs.writeFileSync(fn+'.json',JSON.stringify({pc,dsf,...gt}));k++;}}
 await pg.close();}
 console.log('Zeilen',k);await b.close();})();
