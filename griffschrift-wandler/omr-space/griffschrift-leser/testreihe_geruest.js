window._T = {res: [], log: [], done: false};
(function () {
  const orig = window.omrSend;
  window.omrSend = async function (f, p, o) { const r = await orig(f, p, o); window._T.raw = r; window._T.opt = o; return r; };
})();
const PAGEPIECE = {6:'griffiglandler',7:'griffiglandler',8:'furnvoter',9:'furnvoter',10:'polkaroyal',11:'polkaroyal',12:'zommgschwoasste',13:'zommgschwoasste',14:'eiswind',15:'eiswind',16:'eiswind',17:'isabellile',18:'zehnfuenfer',19:'zehnfuenfer',20:'zehnfuenfer',21:'polarweis',22:'boehmischerhauch',23:'boehmischerhauch',24:'boehmischerhauch',25:'bobblpolka',26:'sibirischegedanken',27:'sibirischegedanken',28:'sibirischegedanken',29:'sibirischegedanken'};
window._TRUTH = {};
async function truth(k) {
  if (!_TRUTH[k]) { const n = JSON.parse(await gunzipB64(PHILIPP[k].p)); const m = {}; n.forEach(([t, pi, r, b, d]) => { const key = Math.round(t * 1000); (m[key] = m[key] || []).push([r, b, d, pi]); }); _TRUTH[k] = Object.keys(m).map(Number).sort((a, b) => a - b).map(t => m[t]); }
  return _TRUTH[k];
}
function chordsFromPins(n) { const m = {}; n.forEach(([t, pi, r, b, d]) => { const key = Math.round(t * 1000); (m[key] = m[key] || []).push([r, b, d, pi]); }); return Object.keys(m).map(Number).sort((a, b) => a - b).map(t => m[t]); }
const S = c => new Set(c.map(x => x[0] + ':' + x[1]));
function sc(a, b) { let i = 0; a.forEach(x => { if (b.has(x)) i++; }); return i ? 2 * i - 0.5 * (a.size + b.size - 2 * i) : -3; }
// semiglobal: A (gelesen) ganz, B (bekannt) nur ein Ausschnitt
function align(A, B) {
  const n = A.length, m = B.length, Sx = [], bt = [];
  for (let i = 0; i <= n; i++) { Sx.push(new Float64Array(m + 1)); bt.push(new Int8Array(m + 1)); Sx[i][0] = -0.5 * i; }
  for (let i = 1; i <= n; i++) for (let j = 1; j <= m; j++) {
    const c = [Sx[i - 1][j - 1] + sc(A[i - 1], B[j - 1]), Sx[i - 1][j] - 0.5, Sx[i][j - 1] - (i === n ? 0 : 0.5)];
    let k = 0; if (c[1] > c[k]) k = 1; if (c[2] > c[k]) k = 2; Sx[i][j] = c[k]; bt[i][j] = k;
  }
  let j = 0, best = -1e9; for (let q = 0; q <= m; q++) if (Sx[n][q] > best) { best = Sx[n][q]; j = q; }
  let i = n; const pr = [];
  while (i > 0 && j > 0) { const k = bt[i][j]; if (k === 0) { pr.push([i - 1, j - 1]); i--; j--; } else if (k === 1) i--; else j--; }
  return pr.reverse();
}
async function waitDone(t0) {
  for (;;) {
    await new Promise(z => setTimeout(z, 1000));
    const s = document.getElementById('scanst').textContent || '';
    if (/^Fertig|Fehler gemeldet|nicht geantwortet|Done\.|did not answer|reported an error/.test(s)) return s;
    if (Date.now() - t0 > 240000) return 'ZEITÜBERSCHREITUNG ' + s;
  }
}
function setSel(id, v) { const e = document.getElementById(id); e.value = v; e.dispatchEvent(new Event('change', {bubbles: true})); }
window.runTests = async function (cases) {
  const inp = document.getElementById('tf'), files = {}; [...inp.files].forEach(f => files[f.name] = f);
  for (const c of cases) {
    const r = {name: c.file, disk: c.disk, stimm: c.stimm};
    try {
      setSel('disk', c.disk); setSel('stimm', c.stimm); setSel('scankind', 'gs');
      if (!document.getElementById('scanok').checked) document.getElementById('scanok').click();
      _T.raw = null; const dt = new DataTransfer(); dt.items.add(files[c.file]); const si = document.getElementById('scanfile'); si.files = dt.files;
      const t0 = Date.now(); si.dispatchEvent(new Event('change', {bubbles: true}));
      r.status = (await waitDone(t0)).slice(0, 160); r.sec = Math.round((Date.now() - t0) / 1000);
      if (_T.raw && /^Fertig|^Done/.test(r.status)) {
        const res = JSON.parse(_T.raw); const n = res.pins.n; r.griffe = n.length; r.r4 = n.filter(x => x[2] === 4).length;
        r.takte = (LASTXML.match(/<measure /g) || []).length; r.sh = res.pins.sh;
        r.layoutGesendet = !!(_T.opt && _T.opt.layout && _T.opt.layout.length > 50);
        r.usepin = document.getElementById('usepin').checked; r.systeme = document.querySelectorAll('#griff svg.gsys').length;
        if (c.piece) {
          const B = await truth(c.piece), A = chordsFromPins(n), pr = align(A.map(S), B.map(S));
          let eq = 0, balg = 0, s4 = 0, g4 = 0, tone = 0, toneN = 0;
          pr.forEach(([i, j]) => {
            const a = S(A[i]), b = S(B[j]); if (a.size === b.size && [...a].every(x => b.has(x))) eq++;
            const bd = B[j][0][2]; if (A[i][0][2] === bd) balg++;
            B[j].forEach(x => { if (x[0] === 4) { s4++; if (a.has(x[0] + ':' + x[1])) g4++; } });
            A[i].forEach(x => { const y = B[j].find(z => z[0] === x[0] && z[1] === x[1] && z[2] === x[2]); if (y) { toneN++; if (y[3] === x[3] + (c.toneOff || 0)) tone++; } });
          });
          Object.assign(r, {akkorde: A.length, gepaart: pr.length, gleich: eq, balg, soll4: s4, gef4: g4, toene: tone, toeneN: toneN});
        }
      }
    } catch (e) { r.status = 'JS-Fehler: ' + e.message; }
    _T.res.push(r);
  }
  _T.done = true;
};
'Testgerüst geladen';
