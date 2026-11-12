"""Documentation routes on the tool server: allowlisted HTML pages and traversal 404.

Run: python3 -m pytest layout_atlas/test_doc_pages.py
"""
from __future__ import annotations

import http.server
import json
import threading
import unittest
import urllib.error
import urllib.request
from pathlib import Path

import doc_pages
import tool_server as ts

DOC_CARDS = [
    'docs/installation.md',
    'docs/data_sources.md',
    'docs/reaction_map_matching.md',
    'docs/methodology.md',
    'docs/layout_evaluation.md',
    'docs/limitations.md',
    'THIRD_PARTY_NOTICES.md',
    'data/README.md',
]


def get(port, path):
    try:
        with urllib.request.urlopen(f'http://127.0.0.1:{port}{path}', timeout=5) as r:
            return r.status, r.headers.get_content_type(), r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.headers.get_content_type() if e.headers else '', e.read()


class DocPagesUnit(unittest.TestCase):
    def test_traversal_rejected(self):
        for rel in ('docs/../LICENSE', 'docs/../../etc/passwd', '../README.md', 'docs/foo/../../LICENSE'):
            self.assertIsNone(doc_pages.normalize(rel))
            self.assertIsNone(doc_pages.resolve(rel))
            self.assertIsNone(doc_pages.handle('/' + rel))

    def test_allowlist_hits(self):
        for rel in DOC_CARDS + ['README.md', 'LICENSE', 'docs/licensing_audit.md', 'docs/ec_extension_recommendation.md']:
            hit = doc_pages.resolve(rel)
            self.assertIsNotNone(hit, rel)
            self.assertTrue(hit[0].is_file(), rel)


class DocPagesServer(unittest.TestCase):
    def setUp(self):
        with ts.lock:
            ts.CLIENTS.clear()
            ts.LIFE.clear()
            ts.LIFE.update(ever=False, empty_since=None, keep=True, pending_exit=False)
        self.httpd = http.server.ThreadingHTTPServer(('127.0.0.1', 0), ts.Handler)
        self.port = self.httpd.server_address[1]
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self.thread.start()

    def tearDown(self):
        with ts.lock:
            ts.LIFE['keep'] = True
        try:
            self.httpd.shutdown()
        except Exception:
            pass
        self.httpd.server_close()

    def test_each_documentation_card_returns_html(self):
        for rel in DOC_CARDS:
            status, ctype, body = get(self.port, '/' + rel)
            self.assertEqual(status, 200, rel)
            self.assertEqual(ctype, 'text/html', rel)
            self.assertIn(b'<html', body.lower())
            self.assertIn(b'class="top"', body)
            self.assertIn(b'Open Meta Maps', body)
            self.assertIn(b'Map Generator', body)
            self.assertIn(b'About / Provenance', body)

    def test_doc_table_css_avoids_midword_breaks(self):
        status, _, body = get(self.port, '/docs/methodology.md')
        self.assertEqual(status, 200)
        html = body.decode()
        self.assertIn('article.doc th,article.doc td{word-break:normal;overflow-wrap:break-word}', html)
        self.assertIn('th:first-child,article.doc td:first-child{width:1%;white-space:nowrap', html)
        self.assertIn('class="tablewrap"', html)
        self.assertIn('<table>', html)
        self.assertIn('EC currentisation', html)

    def test_inline_about_docs_and_image(self):
        for rel in ('docs/licensing_audit.md', 'docs/ec_extension_recommendation.md', 'LICENSE', 'README.md'):
            status, ctype, body = get(self.port, '/' + rel)
            self.assertEqual(status, 200, rel)
            self.assertEqual(ctype, 'text/html', rel)
            self.assertIn(b'class="top"', body)
        status, ctype, body = get(self.port, '/docs/images/workflow.svg')
        self.assertEqual(status, 200)
        self.assertIn('svg', ctype)
        self.assertTrue(body.startswith(b'<svg') or b'<svg' in body[:200])

    def test_path_traversal_returns_404(self):
        for path in (
            '/docs/../LICENSE',
            '/docs/../../etc/passwd',
            '/docs/%2e%2e/LICENSE',
            '/docs/foo/../../README.md',
        ):
            status, ctype, body = get(self.port, path)
            self.assertEqual(status, 404, path)
            if ctype == 'application/json':
                payload = json.loads(body.decode())
                self.assertEqual(payload.get('error', {}).get('code'), 'not_found', path)

    def test_non_allowlisted_repo_path_404(self):
        status, _, body = get(self.port, '/docs/evaluation/summary.csv')
        self.assertEqual(status, 404)
        self.assertIn(b'not_found', body)
