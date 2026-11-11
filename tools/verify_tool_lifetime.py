#!/usr/bin/env python3
"""Headless checks for tool-server lifetime and downloads.

Run: python3 tools/verify_tool_lifetime.py
"""
from __future__ import annotations

import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / 'layout_atlas'
CHROME = shutil.which('google-chrome') or shutil.which('chromium') or shutil.which('chromium-browser')


def wait_health(port, timeout=15, expect_alive=True):
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(f'http://127.0.0.1:{port}/api/health', timeout=1) as r:
                body = json.loads(r.read().decode())
                if expect_alive:
                    return body
        except Exception:
            if not expect_alive:
                return None
        time.sleep(0.2)
    raise SystemExit(f'health wait failed (expect_alive={expect_alive}) on port {port}')


def start_server(port, keep=False):
    cmd = [sys.executable, str(ATLAS / 'tool_server.py'), str(port)]
    if keep:
        cmd.append('--keep-running')
    env = dict(os.environ)
    env['PYTHONPATH'] = str(ATLAS) + os.pathsep + env.get('PYTHONPATH', '')
    proc = subprocess.Popen(cmd, cwd=str(ATLAS), env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    try:
        wait_health(port, timeout=20)
    except SystemExit:
        out = ''
        try:
            proc.kill()
            out = proc.communicate(timeout=2)[0] or ''
        except Exception:
            pass
        raise SystemExit(f'server failed to start on {port}: {out[-500:]}')
    return proc


def stop(proc):
    if proc.poll() is None:
        proc.send_signal(signal.SIGTERM)
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


def post(port, path):
    urllib.request.urlopen(urllib.request.Request(f'http://127.0.0.1:{port}{path}', data=b'', method='POST'), timeout=3)


def chrome_js(port, expression, timeout_ms=20000):
    if not CHROME:
        raise SystemExit('Chrome not found')
    try:
        import websocket  # type: ignore
    except ImportError:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'websocket-client'])
        import websocket  # type: ignore

    td = tempfile.mkdtemp(prefix='opm-chrome-')
    user_data = Path(td) / 'chrome'
    user_data.mkdir()
    debug_port = 9333 + (port % 50)
    chrome = subprocess.Popen([
        CHROME, f'--remote-debugging-port={debug_port}', '--headless=new', '--disable-gpu',
        '--no-first-run', '--no-default-browser-check', f'--user-data-dir={user_data}',
        f'http://127.0.0.1:{port}/tool/index.html',
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        deadline = time.time() + 20
        ws = None
        while time.time() < deadline and not ws:
            try:
                with urllib.request.urlopen(f'http://127.0.0.1:{debug_port}/json/list', timeout=1) as r:
                    for t in json.loads(r.read().decode()):
                        if 'tool/index.html' in t.get('url', ''):
                            ws = t.get('webSocketDebuggerUrl')
                            break
            except Exception:
                time.sleep(0.2)
        if not ws:
            raise SystemExit('no Chrome tab for tool page')
        conn = websocket.create_connection(ws, timeout=10)
        msg_id = 0

        def send(method, params=None):
            nonlocal msg_id
            msg_id += 1
            conn.send(json.dumps({'id': msg_id, 'method': method, 'params': params or {}}))
            while True:
                data = json.loads(conn.recv())
                if data.get('id') == msg_id:
                    return data

        send('Runtime.enable')
        time.sleep(1.0)
        wrapped = (
            '(async function(){ try { const result = await (async function(){ %s })();'
            ' return {ok:true, result}; } catch(e){ return {ok:false, error: String(e && e.stack || e)}; } })()'
        ) % expression
        data = send('Runtime.evaluate', {
            'expression': wrapped,
            'awaitPromise': True,
            'returnByValue': True,
            'timeout': timeout_ms,
        })
        conn.close()
        val = (((data.get('result') or {}).get('result') or {}).get('value'))
        if not val or not val.get('ok'):
            raise SystemExit('Chrome JS failed: ' + json.dumps(data)[:800])
        return val['result']
    finally:
        chrome.kill()
        try:
            chrome.wait(timeout=3)
        except Exception:
            pass
        shutil.rmtree(td, ignore_errors=True)


def main():
    sys.path.insert(0, str(ATLAS))
    import tool_server as ts
    assert ts.STALE >= 180.0, ts.STALE
    assert ts.GRACE <= 30.0, ts.GRACE

    results = []

    # --- keep-running ---
    port = 8769
    proc = start_server(port, keep=True)
    try:
        h = wait_health(port)
        assert h.get('auto_exit') is False, h
        post(port, '/api/ping?c=k1')
        post(port, '/api/bye?c=k1')
        time.sleep(10)
        wait_health(port)
        results.append(('keep-running stays up + auto_exit false', True))
    finally:
        stop(proc)
        time.sleep(0.5)

    # --- API lifetime (no browser clients left behind) ---
    proc = start_server(port, keep=False)
    try:
        post(port, '/api/ping?c=bg')
        time.sleep(50)  # old STALE was 45 s
        wait_health(port)
        results.append(('throttled/hidden silence 50s keeps server (STALE>>45)', True))

        for name in ('map.svg', 'job.json'):
            with urllib.request.urlopen(f'http://127.0.0.1:{port}/tool/mock/{name}') as r:
                assert r.status == 200 and len(r.read()) > 20
        results.append(('SVG and JSON fetches work', True))

        post(port, '/api/ping?c=rel')
        post(port, '/api/bye?c=rel')
        post(port, '/api/bye?c=bg')
        time.sleep(1)
        post(port, '/api/ping?c=rel')
        time.sleep(10)
        wait_health(port)
        results.append(('reload within GRACE keeps server', True))

        post(port, '/api/bye?c=rel')
        wait_health(port, timeout=20, expect_alive=False)
        assert proc.wait(timeout=5) is not None
        results.append(('last tab bye exits after GRACE', True))
    finally:
        stop(proc)
        time.sleep(0.5)

    # --- Chrome page behaviour on a fresh keep-running server ---
    if CHROME:
        proc = start_server(port, keep=True)
        try:
            info = chrome_js(port, r"""
              const beacons = [];
              const orig = navigator.sendBeacon.bind(navigator);
              navigator.sendBeacon = function(url, data){ beacons.push(String(url)); return orig(url, data); };
              let pagehide = 0;
              window.addEventListener('pagehide', function(){ pagehide++; });
              const a = document.createElement('a');
              a.href = '/api/health';
              a.setAttribute('download', 'health.json');
              document.body.appendChild(a);
              a.click();
              await new Promise(r => setTimeout(r, 400));
              document.dispatchEvent(new Event('visibilitychange'));
              window.dispatchEvent(new Event('focus'));
              return {pagehide, beacons, hasVis: !!window.__opmHB};
            """)
            assert info['pagehide'] == 0, info
            assert not any('/api/bye' in b for b in info['beacons']), info
            assert info['hasVis'], info
            results.append(('download click does not send bye / pagehide; heartbeat installed', True))
        finally:
            stop(proc)
    else:
        results.append(('download click (skipped: no Chrome)', True))

    print('verify_tool_lifetime:')
    for name, ok in results:
        print(('PASS' if ok else 'FAIL'), name)
    if not all(ok for _, ok in results):
        raise SystemExit(1)
    print('ALL PASS')


if __name__ == '__main__':
    main()
