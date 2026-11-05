// Trainingsbilder für den Kopf-Erkenner (YOLO): der Wandler setzt Stücke als Griffschrift, je System ein PNG und eine Liste der Köpfe
// mit Rahmen und Form. Die Formen werden zufällig abgewandelt, damit das Modell alle Schreibweisen kennt:
//   0 runder Kopf voll, 1 runder Kopf hohl, 2 Kreuz, 3 Kreuz im Kreis, 4 Dreieck voll, 5 Dreieck hohl,
//   6 vorgesetztes Doppelkreuz (Rosenzopf; der Kopf rechts davon ist dann rund und zählt als Kreuzkopf)
// Für das Modell zählt nur die Form, nicht Reihe oder Ton; Kreuze werden darum zufällig zu Dreiecken oder zu runden Köpfen mit
// vorgesetztem Doppelkreuz umgezeichnet. Seitenaufbau und Fotofehler macht yolo_daten.py.
// Aufruf (im Ordner griffschrift-wandler einen Webserver starten: python3 -m http.server 8798):
//   node yolo_erzeuge_wandlerbilder.js AUSGABEORDNER STARTWERT ANZAHL_STUECKE
const {chromium}=require('playwright');const fs=require('fs');
const [out,seed0,N]=process.argv.slice(2);fs.mkdirSync(out,{recursive:true});
const URL=process.env.WANDLER||'http://127.0.0.1:8798/index.html?lang=de';
(async()=>{const b=await chromium.launch();let k=0;let seed=+seed0;
 const rnd=()=>{seed=(seed*9301+49297)%233280;return seed/233280;};
 const pg0=await (await b.newContext()).newPage();await pg0.goto(URL);await pg0.waitForTimeout(2500);
 const VM=await pg0.evaluate(()=>[...document.querySelectorAll('#demo option')].map(o=>o.value).filter(v=>v.startsWith('vm:')||v.startsWith('mid:')));
 await pg0.close();
 for(const dsf of [1.5,2,2.5]){
 const pg=await (await b.newContext({locale:'de-DE',viewport:{width:1100+Math.round(rnd()*300),height:1400},deviceScaleFactor:dsf,colorScheme:'light'})).newPage();
 pg.on('dialog',d=>d.dismiss());
 await pg.goto(URL);await pg.waitForTimeout(2500);
 for(let ii=0;ii<+N;ii++){const pc=VM[Math.floor(rnd()*VM.length)];
  const vo=['1','2','3','4'][Math.floor(rnd()*4)];
  try{
   await pg.evaluate(async([pc,vo])=>{const set=(id,v)=>{const s=document.getElementById(id);s.value=v;s.dispatchEvent(new Event('change',{bubbles:true}));};
     set('voices',vo);set('view','griff');const [kind,key]=pc.split(':');if(kind==='vm')await loadVM(key);else{const ds=document.getElementById('demo');ds.value=pc;await loadDemo(pc);}},[pc,vo]);
  }catch(e){continue;}
  await pg.waitForTimeout(1500);
  const n=await pg.evaluate(()=>document.querySelectorAll('#griff svg.gsys').length);
  // je Stück eine Schreibweise: 0 Fuchs (Kreuze), 1 Margetich (4. Reihe Dreieck), 2 Rosenzopf (vorgesetzte Doppelkreuze), 3 gemischt
  const art=Math.floor(rnd()*4);
  for(let i=0;i<n;i++){
   const gt=await pg.evaluate(([i,sd,art])=>{let s0=sd;const rnd=()=>{s0=(s0*9301+49297)%233280;return s0/233280;};
    const NS='http://www.w3.org/2000/svg';
    const s=document.querySelectorAll('#griff svg.gsys')[i];
    const sl=[...s.querySelectorAll('line.sl')].filter(l=>l.getAttribute('y1')===l.getAttribute('y2')).map(l=>+l.getAttribute('y1'));
    const ys=[...new Set(sl)].sort((a,b)=>a-b);if(ys.length<5)return null;const d=(ys[4]-ys[0])/4;
    const gray=Math.round(20+150*rnd());s.querySelectorAll('line.sl').forEach(l=>{l.style.stroke=`rgb(${gray},${gray},${gray})`;l.style.strokeWidth=0.6+0.9*rnd();});
    s.querySelectorAll('text.mn').forEach(t=>t.style.display=rnd()<0.5?'none':'');s.querySelectorAll('rect.brk,rect.msel,rect.bhit,rect.chit,circle.hh').forEach(t=>t.style.display='none');   // Klickflächen nie zeigen
    document.querySelectorAll('body *').forEach(el=>{if(getComputedStyle(el).position==='fixed')el.style.display='none';});
    s.querySelectorAll('[style]').forEach(el=>{const st=el.getAttribute('style');if(/var\(--(extra|muted|accent|zug|druck)\)/.test(st))el.setAttribute('style',st.replace(/var\(--(extra|muted|accent|zug|druck)\)/g,'#000'));});
    s.querySelectorAll('.cur,.selx,.forced').forEach(el=>el.classList.remove('cur','selx','forced'));
    const es=0.85+0.3*rnd(),xs=0.8+0.35*rnd(),xsw=1.2+1.4*rnd(),tw=(1.1+0.4*rnd())*d,th=(0.9+0.3*rnd())*d,tsw=0.8+1.0*rnd();
    const items=[];   // [Element, Klasse]
    s.querySelectorAll('g.gn').forEach(g=>{
      [...g.querySelectorAll('ellipse')].forEach(el=>{el.setAttribute('rx',5.4*es);el.setAttribute('ry',3.9*es);const st=el.getAttribute('style')||'';items.push([el,/fill:none/.test(st)?1:0]);});
      [...g.querySelectorAll('path')].forEach(p=>{const m=(p.getAttribute('d')||'').match(/^M([-\d.]+) ([-\d.]+) L([-\d.]+) ([-\d.]+) M/);if(!m)return;
        const cx=(+m[1]+ +m[3])/2,cy=(+m[2]+ +m[4])/2,hollow=+(p.getAttribute('stroke-width')||2)<2;
        const circ=[...g.querySelectorAll('circle:not(.hh)')].find(c=>Math.abs(+c.getAttribute('cx')-cx)<0.6&&Math.abs(+c.getAttribute('cy')-cy)<0.6);
        const wahl=art===1?(rnd()<0.6?'t':'x'):art===2?(hollow?'x':(rnd()<0.75?'p':'x')):art===3?['x','t','p'][Math.floor(rnd()*3)]:'x';
        if(wahl==='t'){   // Dreieck, Spitze oben, hohl ab der Halben
          if(circ)circ.remove();const poly=document.createElementNS(NS,'polygon');
          poly.setAttribute('points',`${cx-tw/2},${cy+th/2} ${cx+tw/2},${cy+th/2} ${cx},${cy-th/2}`);
          poly.setAttribute('style',hollow?`fill:none;stroke:#000;stroke-width:${tsw}`:'fill:#000;stroke:#000;stroke-width:0.5');
          p.replaceWith(poly);items.push([poly,hollow?5:4]);}
        else if(wahl==='p'){   // runder Kopf mit vorgesetztem Doppelkreuz (Rosenzopf)
          if(circ)circ.remove();const e=document.createElementNS(NS,'ellipse');e.setAttribute('cx',cx);e.setAttribute('cy',cy);e.setAttribute('rx',5.4*es);e.setAttribute('ry',3.9*es);
          e.setAttribute('transform',`rotate(-20 ${cx} ${cy})`);e.setAttribute('style','fill:#000');p.replaceWith(e);items.push([e,0]);
          const gx=cx-(1.25+0.35*rnd())*d,r=0.32*d,gp=document.createElementNS(NS,'g');
          const ln=(x1,y1,x2,y2,w)=>{const l=document.createElementNS(NS,'line');l.setAttribute('x1',x1);l.setAttribute('y1',y1);l.setAttribute('x2',x2);l.setAttribute('y2',y2);l.setAttribute('style',`stroke:#000;stroke-width:${w}`);gp.appendChild(l);};
          const w0=0.18*d+0.1*d*rnd();ln(gx-r,cy-r,gx+r,cy+r,w0);ln(gx-r,cy+r,gx+r,cy-r,w0);
          for(const [sx,sy] of [[-1,-1],[1,-1],[-1,1],[1,1]]){const q=document.createElementNS(NS,'rect');const a=(0.24+0.08*rnd())*d;q.setAttribute('x',gx+sx*r-a/2);q.setAttribute('y',cy+sy*r-a/2);q.setAttribute('width',a);q.setAttribute('height',a);q.setAttribute('style','fill:#000');gp.appendChild(q);}
          g.appendChild(gp);items.push([gp,6]);}
        else{const rr=4.5*xs;p.setAttribute('d',`M${cx-rr} ${cy-rr} L${cx+rr} ${cy+rr} M${cx-rr} ${cy+rr} L${cx+rr} ${cy-rr}`);p.setAttribute('stroke-width',hollow?1.3:xsw);
          if(hollow&&!circ){const c=document.createElementNS(NS,'circle');c.setAttribute('cx',cx);c.setAttribute('cy',cy);c.setAttribute('r',0.62*d);c.setAttribute('style','fill:none;stroke:#000;stroke-width:1');g.appendChild(c);}
          items.push([p,hollow?3:2]);}});});
    const R=s.getBoundingClientRect();
    const boxes=items.map(([el,c])=>{const r=el.getBoundingClientRect();return [c,r.left-R.left,r.top-R.top,r.width,r.height];}).filter(q=>q[3]>1&&q[4]>1);
    return {boxes,w:R.width,h:R.height,staff:ys};},[i,Math.floor(rnd()*100000),art]);
   if(!gt||!gt.boxes.length)continue;
   const el=(await pg.$$('#griff svg.gsys'))[i];const fn=`${out}/y${seed0}_${String(k).padStart(5,'0')}`;await el.screenshot({path:fn+'.png'});
   fs.writeFileSync(fn+'.json',JSON.stringify({pc,dsf,art,...gt}));k++;}}
 await pg.close();}
 console.log('Systeme',k);await b.close();})();
