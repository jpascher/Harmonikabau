// Wandler-Bilder mit Bassbuchstaben als Lerndaten: je System ein PNG und die Buchstaben mit Rechteck (Bildpixel)
// Aufruf: node gen_bass.js OUTDIR SEED N   (N Stücke aus midi/vm)
const {chromium}=require('playwright');const fs=require('fs');
const [out,seed0,N]=process.argv.slice(2);fs.mkdirSync(out,{recursive:true});
const vm=fs.readdirSync('/tmp/claude-0/srv/midi/vm').filter(f=>f.endsWith('.mid')).map(f=>f.slice(0,-4));
(async()=>{const b=await chromium.launch();let seed=+seed0;const rnd=()=>{seed=(seed*9301+49297)%233280;return seed/233280;};let k=0;
 for(const dsf of [1,1.5,2]){
  const pg=await (await b.newContext({locale:'de-DE',viewport:{width:1000+Math.round(rnd()*400),height:1400},deviceScaleFactor:dsf,colorScheme:'light'})).newPage();
  pg.on('dialog',d=>d.dismiss());
  await pg.goto('http://127.0.0.1:8798/index.html?lang=de');await pg.waitForTimeout(2500);
  const styles=await pg.evaluate(()=>[...document.getElementById('gstyle').options].map(o=>o.value));
  for(let i=0;i<+N;i++){const key=vm[Math.floor(rnd()*vm.length)];const sty=styles[Math.floor(rnd()*styles.length)];
   try{await pg.evaluate(async([key,sty])=>{const set=(id,v)=>{const s=document.getElementById(id);if(!s)return;s.value=v;s.dispatchEvent(new Event('change',{bubbles:true}));};
     set('gstyle',sty);set('view','griff');await loadVM(key);},[key,sty]);}catch(e){continue;}
   await pg.waitForTimeout(1200);
   const n=await pg.evaluate(()=>document.querySelectorAll('#griff svg.gsys').length);
   for(let j=0;j<n;j++){
    const gt=await pg.evaluate(([j,dsf])=>{const s=document.querySelectorAll('#griff svg.gsys')[j];const R=s.getBoundingClientRect();
     document.querySelectorAll('body *').forEach(el=>{if(getComputedStyle(el).position==='fixed')el.style.display='none';});
     s.querySelectorAll('.cur,.selx,.forced,.miss').forEach(el=>el.classList.remove('cur','selx','forced','miss'));
     s.querySelectorAll('[style]').forEach(el=>{const st=el.getAttribute('style');if(/var\(--/.test(st))el.setAttribute('style',st.replace(/var\(--[a-z-]+\)/g,'#000'));});
     s.querySelectorAll('text,line,path,ellipse,polygon,circle,rect').forEach(el=>{const c=getComputedStyle(el);if(c.fill&&c.fill!=='none'&&!/rgba\(0, 0, 0, 0\)/.test(c.fill))el.style.fill='#000';if(c.stroke&&c.stroke!=='none')el.style.stroke='#000';});
     s.querySelectorAll('rect.brk,rect.msel,rect.bhit,rect.chit,rect.selbg').forEach(el=>el.style.display='none');
     const sl=[...s.querySelectorAll('line.sl')].filter(l=>l.getAttribute('y1')===l.getAttribute('y2')).map(l=>+l.getAttribute('y1'));
     const L=[...s.querySelectorAll('text.bl')].map(t=>{const r=t.getBoundingClientRect();return {t:t.textContent,x0:(r.left-R.left)*dsf,y0:(r.top-R.top)*dsf,x1:(r.right-R.left)*dsf,y1:(r.bottom-R.top)*dsf};}).filter(x=>x.t.trim());
     return {L,staff:[...new Set(sl)].sort((a,b)=>a-b).map(y=>y*dsf)};},[j,dsf]);
    if(!gt.L.length||gt.staff.length<5)continue;
    const el=(await pg.$$('#griff svg.gsys'))[j];const fn=`${out}/w${seed0}_${String(k).padStart(4,'0')}`;await el.screenshot({path:fn+'.png'});
    fs.writeFileSync(fn+'.json',JSON.stringify({key,sty,dsf,...gt}));k++;}}
  await pg.close();}
 console.log('Systeme',k);await b.close();})();
