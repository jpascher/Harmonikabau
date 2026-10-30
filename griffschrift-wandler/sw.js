// Griffschrift-Wandler: Programm offline verfügbar halten (Seite: zuerst Netz, sonst Zwischenspeicher; MIDI-Dateien: Zwischenspeicher)
const C="gw-v1";
self.addEventListener("install",e=>{e.waitUntil(caches.open(C).then(c=>c.addAll(["./","index.html","manifest.webmanifest","icon-192.png","icon-512.png"])).then(()=>self.skipWaiting()));});
self.addEventListener("activate",e=>{e.waitUntil(caches.keys().then(ks=>Promise.all(ks.filter(k=>k!==C).map(k=>caches.delete(k)))).then(()=>self.clients.claim()));});
self.addEventListener("fetch",e=>{const u=new URL(e.request.url);if(e.request.method!=="GET"||u.origin!==location.origin)return;
  if(/\/midi\//.test(u.pathname)){e.respondWith(caches.open(C).then(c=>c.match(e.request).then(r=>r||fetch(e.request).then(n=>{if(n.ok)c.put(e.request,n.clone());return n;}))));return;}
  e.respondWith(fetch(e.request).then(n=>{if(n.ok){const cl=n.clone();caches.open(C).then(c=>c.put(e.request,cl));}return n;}).catch(()=>caches.match(e.request).then(r=>r||caches.match("index.html"))));});
