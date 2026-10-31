// Griffschrift-Wandler: Programm offline verfügbar halten (Seite: zuerst Netz, sonst Zwischenspeicher; MIDI- und Klang-Dateien: Zwischenspeicher)
const C="gw-v2";   // Seite immer beim Server nachfragen (no-cache), damit ein Push sofort sichtbar wird
self.addEventListener("install",e=>{e.waitUntil(caches.open(C).then(c=>c.addAll(["./","index.html","manifest.webmanifest","icon-192.png","icon-512.png"])).then(()=>self.skipWaiting()));});
self.addEventListener("activate",e=>{e.waitUntil(caches.keys().then(ks=>Promise.all(ks.filter(k=>k!==C).map(k=>caches.delete(k)))).then(()=>self.clients.claim()));});
self.addEventListener("fetch",e=>{const u=new URL(e.request.url);if(e.request.method!=="GET"||u.origin!==location.origin)return;
  if(/\/(midi|klang)\//.test(u.pathname)){e.respondWith(caches.open(C).then(c=>c.match(e.request).then(r=>r||fetch(e.request).then(n=>{if(n.ok)c.put(e.request,n.clone());return n;}))));return;}
  // zuerst Netz, aber höchstens 4 s warten (schlechter Handyempfang), dann den gespeicherten Stand zeigen
  const net=fetch(u.href,{cache:"no-cache"}).then(n=>{if(n.ok){const cl=n.clone();caches.open(C).then(c=>c.put(e.request,cl));}return n;});
  const old=()=>caches.match(e.request).then(r=>r||caches.match("index.html"));
  e.respondWith(new Promise(ok=>{let done=false;const fin=r=>{if(!done&&r){done=true;ok(r);}};
    const t=setTimeout(()=>old().then(fin),4000);
    net.then(r=>{clearTimeout(t);fin(r);}).catch(()=>{clearTimeout(t);old().then(r=>fin(r||Response.error()));});}));});
