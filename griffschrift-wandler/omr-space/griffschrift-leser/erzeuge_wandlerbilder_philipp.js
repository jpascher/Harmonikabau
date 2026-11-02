// Aufruf: node erzeuge_wandlerbilder_philipp.js griffiglandler furnvoter …  (Bilder nach /tmp/philipp/wand)
// Wandler-Griffschrift (Fuchs-Schreibweise: 3. und 4. Reihe als Kreuze) von Philipps Stücken als Bilder zeichnen
const {chromium} = require('playwright'); const fs = require('fs');
(async () => {
  const b = await chromium.launch();
  const pg = await (await b.newContext({locale: 'de-DE', viewport: {width: 1100, height: 1400}, deviceScaleFactor: 2, colorScheme: 'light'})).newPage();
  pg.on('dialog', d => d.dismiss());
  await pg.goto('file://'+require('path').resolve(__dirname,'../../index.html')+'?lang=de'); await pg.waitForTimeout(2500);
  const out = {};
  for (const k of process.argv.slice(2)) {
    await pg.evaluate(async k => { const s = document.getElementById('demo'); s.value = 'phx:' + k; s.dispatchEvent(new Event('change')); }, k);
    await pg.waitForTimeout(2500);
    const n = await pg.evaluate(() => document.querySelectorAll('#griff svg.gsys').length);
    out[k] = n;
    for (let i = 0; i < n; i++) {
      const el = (await pg.$$('#griff svg.gsys'))[i];
      await el.screenshot({path: `/tmp/philipp/wand/${k}_${i}.png`});
    }
  }
  console.log(JSON.stringify(out)); await b.close();
})();
