#!/usr/bin/env node
/**
 * Headless Chrome check: edit mode on a finished tool job.
 * Usage: node test_map_edit_chrome.js <baseUrl> <jobId> [screenshotPath]
 */
const fs = require('fs');
const path = require('path');
const pwPath = [
  path.join(__dirname, 'runtime/node_modules/playwright-core'),
  '/tmp/node_modules/playwright-core',
].find((p) => fs.existsSync(p));
const { chromium } = require(pwPath);

async function main() {
  const base = process.argv[2] || 'http://127.0.0.1:8765';
  const jobId = process.argv[3];
  const shot = process.argv[4] || '';
  if (!jobId) {
    console.error('usage: node test_map_edit_chrome.js <baseUrl> <jobId> [screenshot]');
    process.exit(2);
  }
  const chrome =
    process.env.CHROME_PATH ||
    ['/usr/local/bin/google-chrome', '/usr/bin/google-chrome-stable', '/usr/bin/chromium'].find((p) =>
      fs.existsSync(p)
    );
  const browser = await chromium.launch({
    executablePath: chrome,
    headless: true,
    args: ['--no-sandbox', '--disable-gpu'],
  });
  const page = await browser.newPage({ viewport: { width: 1280, height: 900 } });
  const url = `${base}/tool/index.html?job=${jobId}`;
  page.setDefaultTimeout(120000);
  await page.goto(url, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.OPM_SVG && window.OPM_PZ && window.OPM_EDIT, null, { timeout: 120000 });
  await page.waitForFunction(() => window.OPM_EDIT && window.OPM_EDIT.isReady && window.OPM_EDIT.isReady() && document.getElementById('pedit'), null, {
    timeout: 120000,
  });
  await page.waitForTimeout(500);

  const before = await page.evaluate(async (id) => {
    const r = await fetch('/api/jobs/' + id + '/files/map.json');
    const j = await r.json();
    const nid = Object.keys(j[1].nodes).find((k) => j[1].nodes[k].node_type === 'metabolite');
    const n = j[1].nodes[nid];
    let original = null;
    try {
      const o = await fetch('/api/jobs/' + id + '/files/map.original.json');
      if (o.ok) original = await o.json();
    } catch (e) {}
    return {
      nid,
      x: n.x,
      y: n.y,
      label_x: n.label_x,
      label_y: n.label_y,
      originalX: original ? original[1].nodes[nid].x : null,
      originalY: original ? original[1].nodes[nid].y : null,
      hadOriginal: !!original,
    };
  }, jobId);

  await page.click('#pedit');
  await page.waitForFunction(() => window.OPM_EDIT.isActive());

  // Drag one metabolite node by ~40 px in screen space
  const dragInfo = await page.evaluate((nid) => {
    const svg = window.OPM_SVG;
    const c = svg.querySelector('circle[data-n="' + CSS.escape(nid) + '"]');
    if (!c) return null;
    const r = c.getBoundingClientRect();
    return { x: r.left + r.width / 2, y: r.top + r.height / 2, cx: +c.getAttribute('cx'), cy: +c.getAttribute('cy') };
  }, before.nid);
  if (!dragInfo) throw new Error('metabolite circle not found: ' + before.nid);

  await page.mouse.move(dragInfo.x, dragInfo.y);
  await page.mouse.down();
  await page.mouse.move(dragInfo.x + 50, dragInfo.y + 30, { steps: 8 });
  await page.mouse.up();
  await page.waitForTimeout(200);

  const afterDrag = await page.evaluate((nid) => {
    const c = window.OPM_SVG.querySelector('circle[data-n="' + CSS.escape(nid) + '"]');
    return { cx: +c.getAttribute('cx'), cy: +c.getAttribute('cy'), dirty: window.OPM_EDIT.isDirty() };
  }, before.nid);
  if (!(afterDrag.cx !== dragInfo.cx || afterDrag.cy !== dragInfo.cy)) throw new Error('drag did not move node');
  if (!afterDrag.dirty) throw new Error('expected dirty after drag');

  // Undo
  await page.click('#pendo');
  await page.waitForTimeout(150);
  const afterUndo = await page.evaluate((nid) => {
    const c = window.OPM_SVG.querySelector('circle[data-n="' + CSS.escape(nid) + '"]');
    return { cx: +c.getAttribute('cx'), cy: +c.getAttribute('cy') };
  }, before.nid);
  if (Math.abs(afterUndo.cx - dragInfo.cx) > 0.2 || Math.abs(afterUndo.cy - dragInfo.cy) > 0.2) {
    throw new Error('undo did not restore node: ' + JSON.stringify({ afterUndo, dragInfo }));
  }

  // Re-drag, then box-select move another region
  await page.mouse.move(dragInfo.x, dragInfo.y);
  await page.mouse.down();
  await page.mouse.move(dragInfo.x + 40, dragInfo.y + 20, { steps: 6 });
  await page.mouse.up();
  await page.waitForTimeout(150);

  // Rubber-band around several on-screen metabolite circles
  const box = await page.evaluate(() => {
    const circles = [...window.OPM_SVG.querySelectorAll('circle[data-n]')].slice(0, 8);
    const rs = circles.map((c) => c.getBoundingClientRect());
    const left = Math.min(...rs.map((r) => r.left)) - 8;
    const top = Math.min(...rs.map((r) => r.top)) - 8;
    const right = Math.max(...rs.map((r) => r.right)) + 8;
    const bottom = Math.max(...rs.map((r) => r.bottom)) + 8;
    return { x: left, y: top, x2: right, y2: bottom };
  });
  await page.mouse.move(box.x, box.y);
  await page.mouse.down();
  await page.mouse.move(box.x2, box.y2, { steps: 6 });
  await page.mouse.up();
  await page.waitForTimeout(100);
  const selCount = await page.evaluate(() => window.OPM_SVG.querySelectorAll('.sel').length);
  if (selCount < 2) throw new Error('box-select selected too few elements: ' + selCount);
  const mid = await page.evaluate(() => {
    const el = window.OPM_SVG.querySelector('.sel');
    const r = el.getBoundingClientRect();
    return { x: r.left + r.width / 2, y: r.top + r.height / 2 };
  });
  await page.mouse.move(mid.x, mid.y);
  await page.mouse.down();
  await page.mouse.move(mid.x + 25, mid.y + 15, { steps: 5 });
  await page.mouse.up();
  await page.waitForTimeout(100);

  if (shot) {
    await page.screenshot({ path: shot, fullPage: false });
  }

  // Save
  await page.click('#psave');
  await page.waitForFunction(() => !window.OPM_EDIT.isDirty(), null, { timeout: 120000 });
  await page.waitForTimeout(500);

  const saved = await page.evaluate(async (id) => {
    const r = await fetch('/api/jobs/' + id + '/files/map.json', { cache: 'no-store' });
    return r.json();
  }, jobId);
  const sn = saved[1].nodes[before.nid];
  if (Math.abs(sn.x - before.x) < 0.5 && Math.abs(sn.y - before.y) < 0.5) {
    throw new Error('saved map.json positions unchanged');
  }

  // Reload page and confirm persistence
  await page.goto(url, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.OPM_SVG && window.OPM_EDIT, null, { timeout: 120000 });
  await page.waitForTimeout(600);
  const reloaded = await page.evaluate((nid) => {
    const c = window.OPM_SVG.querySelector('circle[data-n="' + CSS.escape(nid) + '"]');
    return c ? { cx: +c.getAttribute('cx'), cy: +c.getAttribute('cy') } : null;
  }, before.nid);
  if (!reloaded) throw new Error('node missing after reload');
  if (Math.abs(reloaded.cx - sn.x) > 1 || Math.abs(reloaded.cy - sn.y) > 1) {
    throw new Error('reload SVG mismatch: ' + JSON.stringify({ reloaded, sn }));
  }

  // Download SVG/JSON include edits
  const dlJson = await page.evaluate(async (id) => {
    const r = await fetch('/api/jobs/' + id + '/files/map.json');
    return r.json();
  }, jobId);
  const dlSvg = await page.evaluate(async (id) => {
    const r = await fetch('/api/jobs/' + id + '/files/map.svg');
    return r.text();
  }, jobId);
  if (Math.abs(dlJson[1].nodes[before.nid].x - sn.x) > 0.01) throw new Error('download JSON mismatch');
  if (!dlSvg.includes('cx="' + Number(sn.x).toFixed(1) + '"') && !dlSvg.includes('cx="' + sn.x)) {
    // tolerate formatting
    const m = dlSvg.match(new RegExp('data-n="' + before.nid + '"[^>]*cx="([0-9.]+)"'));
    if (!m || Math.abs(+m[1] - sn.x) > 1) throw new Error('download SVG missing edited cx');
  }

  // original kept (frozen at first save; must not change during this save)
  const orig = await page.evaluate(async (id) => {
    const r = await fetch('/api/jobs/' + id + '/files/map.original.json');
    if (!r.ok) return null;
    return r.json();
  }, jobId);
  if (!orig) throw new Error('map.original.json missing');
  const ox = orig[1].nodes[before.nid].x;
  const oy = orig[1].nodes[before.nid].y;
  if (before.hadOriginal) {
    if (Math.abs(ox - before.originalX) > 0.01 || Math.abs(oy - before.originalY) > 0.01) {
      throw new Error('original was overwritten');
    }
  } else if (Math.abs(ox - before.x) > 0.01 || Math.abs(oy - before.y) > 0.01) {
    throw new Error('first original snapshot does not match pre-edit map');
  }

  console.log(JSON.stringify({ ok: true, jobId, node: before.nid, before: { x: before.x, y: before.y }, after: { x: sn.x, y: sn.y }, selCount }, null, 2));
  await browser.close();
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
