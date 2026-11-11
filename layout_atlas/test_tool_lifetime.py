"""Tool-server lifetime: bye/GRACE primary exit, multi-minute STALE fallback, keep-running.

Run: python3 layout_atlas/test_tool_lifetime.py
"""
from __future__ import annotations

import http.server
import json
import os
import sys
import threading
import time
import unittest
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tool_server as ts  # noqa: E402

TOOL_HTML = Path(__file__).resolve().parent / 'results' / 'tool' / 'index.html'


def reset_life(keep=False):
    with ts.lock:
        ts.CLIENTS.clear()
        ts.LIFE.clear()
        ts.LIFE.update(ever=False, empty_since=None, keep=keep, pending_exit=False)
    ts.WAKE.set()


def post(port, path):
    req = urllib.request.Request(f'http://127.0.0.1:{port}{path}', data=b'', method='POST')
    with urllib.request.urlopen(req, timeout=3) as r:
        return r.status, json.loads(r.read().decode())


def get(port, path):
    with urllib.request.urlopen(f'http://127.0.0.1:{port}{path}', timeout=3) as r:
        return r.status, r.read()


class LifetimeConstants(unittest.TestCase):
    def test_stale_tolerates_background_throttling(self):
        # Chrome intensive throttling ~1 ping/min; several minutes of silence before fallback close.
        self.assertGreaterEqual(ts.STALE, 180.0)
        self.assertLess(ts.GRACE, 30.0)

    def test_heartbeat_pings_on_visibility_and_focus(self):
        hb = ts.Handler.HEARTBEAT.decode()
        self.assertIn("/api/ping?c=", hb)
        self.assertIn("/api/bye?c=", hb)
        self.assertIn("visibilitychange", hb)
        self.assertIn("visibilityState==='visible'", hb)
        self.assertIn("addEventListener('focus'", hb)
        self.assertIn("pagehide", hb)

    def test_tool_page_safe_download_and_visibility_ping(self):
        html = TOOL_HTML.read_text(encoding='utf-8')
        self.assertIn('function safeDownload', html)
        self.assertIn('data-dl=', html)
        self.assertIn('Cannot download: the tool server is not reachable', html)
        self.assertNotIn('无法下载：工具服务器不可用', html)
        self.assertNotIn('id="lang"', html)
        self.assertNotIn('data-zh', html)
        self.assertIn("visibilitychange", html)
        self.assertIn("visibilityState==='visible'", html)
        # bye only on pagehide/unload — not on download clicks
        self.assertIn("addEventListener('pagehide',bye)", html)
        self.assertNotIn('download].forEach', html)  # sanity: keep lean helpers


class IdleWatch(unittest.TestCase):
    def setUp(self):
        reset_life()
        self._grace, self._stale = ts.GRACE, ts.STALE
        ts.GRACE, ts.STALE = 0.4, 1.2
        self.httpd = http.server.ThreadingHTTPServer(('127.0.0.1', 0), ts.Handler)
        self.port = self.httpd.server_address[1]
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self.thread.start()
        self.watch = threading.Thread(target=ts.idle_watch, args=(self.httpd,), daemon=True)
        self.watch.start()

    def tearDown(self):
        # Park any leftover idle_watch threads so they do not observe the next test's CLIENTS.
        with ts.lock:
            ts.LIFE['keep'] = True
        ts.WAKE.set()
        try:
            self.httpd.shutdown()
        except Exception:
            pass
        self.httpd.server_close()
        ts.GRACE, ts.STALE = self._grace, self._stale
        reset_life()

    def alive(self):
        try:
            get(self.port, '/api/health')
            return True
        except Exception:
            return False

    def test_bye_then_grace_exits(self):
        post(self.port, '/api/ping?c=tab1')
        self.assertTrue(self.alive())
        post(self.port, '/api/bye?c=tab1')
        deadline = time.time() + 3.0
        while time.time() < deadline and self.alive():
            time.sleep(0.1)
        self.assertFalse(self.alive())
        self.assertFalse(self.watch.is_alive())

    def test_reload_within_grace_survives(self):
        post(self.port, '/api/ping?c=tab1')
        post(self.port, '/api/bye?c=tab1')
        time.sleep(0.15)
        post(self.port, '/api/ping?c=tab1')  # reload reconnects before GRACE
        time.sleep(ts.GRACE + 0.5)
        self.assertTrue(self.alive())
        self.assertTrue(self.watch.is_alive())

    def test_stale_fallback_exits_after_silence(self):
        post(self.port, '/api/ping?c=tab1')
        # no bye: wait past STALE + GRACE
        deadline = time.time() + ts.STALE + ts.GRACE + 2.0
        while time.time() < deadline and self.alive():
            time.sleep(0.15)
        self.assertFalse(self.alive())

    def test_keep_running_never_exits(self):
        reset_life(keep=True)
        post(self.port, '/api/ping?c=tab1')
        post(self.port, '/api/bye?c=tab1')
        time.sleep(ts.GRACE + 0.8)
        status, body = get(self.port, '/api/health')
        self.assertEqual(status, 200)
        health = json.loads(body.decode())
        self.assertFalse(health.get('auto_exit'))
        self.assertTrue(self.alive())


class ClientHelpers(unittest.TestCase):
    def setUp(self):
        reset_life()

    def tearDown(self):
        reset_life()

    def test_touch_clears_empty_since(self):
        ts.client_touch('a')
        with ts.lock:
            ts.CLIENTS.pop('a', None)
            ts.LIFE['empty_since'] = time.time() - 100
        ts.client_touch('a')
        with ts.lock:
            self.assertIsNone(ts.LIFE['empty_since'])
            self.assertFalse(ts.LIFE['pending_exit'])
            self.assertIn('a', ts.CLIENTS)

    def test_bye_marks_pending_when_last(self):
        ts.client_touch('a')
        ts.client_bye('a')
        with ts.lock:
            self.assertEqual(ts.CLIENTS, {})
            self.assertTrue(ts.LIFE['pending_exit'])


if __name__ == '__main__':
    # Avoid starting the real tool environment check path; Handler only needs module state.
    os.chdir(Path(__file__).resolve().parent)
    unittest.main(verbosity=2)
