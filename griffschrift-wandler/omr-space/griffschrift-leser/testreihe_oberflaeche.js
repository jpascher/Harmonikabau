// Testreihe über die Wandler-Oberfläche: Wandler unter http://localhost:8000 (OMR_URL auf http://localhost:8765 gestellt),
// Erkennungsserver aus omr-space/vps lokal auf Port 8765; Testdateien und cases.json in /tmp/testreihe.
// Aufruf: node testreihe_oberflaeche.js A de
const {chromium} = require('playwright'); const fs = require('fs'); const path = require('path');
(async () => {
  const b = await chromium.launch();
  const pg = await (await b.newContext({locale: 'de-DE', viewport: {width: 1300, height: 1000}})).newPage();
  const errs = []; pg.on('pageerror', e => errs.push(e.message)); pg.on('dialog', d => { errs.push('Dialog: ' + d.message()); d.dismiss(); });
  await pg.goto('http://localhost:8000/griffschrift-wandler/index.html?lang='+(process.argv[3]||'de')); await pg.waitForTimeout(2500);
  await pg.evaluate(() => { const i = document.createElement('input'); i.type = 'file'; i.multiple = true; i.id = 'tf'; document.body.appendChild(i); });
  const dir = '/tmp/testreihe'; const files = fs.readdirSync(dir).filter(f => /\.(png|jpg|pdf)$/.test(f));
  await pg.setInputFiles('#tf', files.map(f => path.join(dir, f)));
  await pg.evaluate(fs.readFileSync(__dirname + '/testreihe_geruest.js', 'utf8'));
  const cases = JSON.parse(fs.readFileSync(dir + '/cases.json', 'utf8'));
  await pg.evaluate(c => { runTests(c); }, cases);
  let last = 0;
  for (;;) {
    await pg.waitForTimeout(3000);
    const st = await pg.evaluate(() => ({n: _T.res.length, done: _T.done}));
    if (st.n !== last) { last = st.n; const r = await pg.evaluate(() => _T.res[_T.res.length - 1]); console.log(JSON.stringify(r)); }
    if (st.done) break;
  }
  const res = await pg.evaluate(() => _T.res);
  fs.writeFileSync(dir + '/ergebnis_' + (process.argv[2] || 'lauf') + '.json', JSON.stringify(res, null, 1));
  // Abspielen nach dem letzten Scan: darf keinen Fehler werfen
  await pg.evaluate(() => { try { document.querySelector('#play, button[id*=play]') && document.querySelector('#play, button[id*=play]').click(); } catch (e) {} });
  await pg.waitForTimeout(1500);
  console.log('Seitenfehler:', JSON.stringify(errs));
  await b.close();
})();
