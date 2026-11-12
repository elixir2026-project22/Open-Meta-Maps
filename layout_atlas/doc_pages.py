"""Allowlisted repository documentation for the local tool server.

Serves docs/**/*.md, README.md, THIRD_PARTY_NOTICES.md, LICENSE, data/README.md
as HTML (site chrome + top nav), and docs/images/* as static files. Relative
markdown links/images that stay inside the allowlist are rewritten to site-root
paths so they work from any doc URL.
"""
from __future__ import annotations

import html
import mimetypes
import re
from pathlib import Path
from urllib.parse import unquote, urlparse

from core import ROOT
import ui

ROOT_DOCS = frozenset({'README.md', 'THIRD_PARTY_NOTICES.md', 'LICENSE', 'data/README.md'})
_H1 = re.compile(r'^#\s+(.+?)\s*$', re.M)
_ATTR = re.compile(r'\b(href|src)=(["\'])([^"\']*)\2', re.I)
# Shared chrome sets overflow-wrap:anywhere on main>article (min-content → one glyph).
# Doc pages restore normal wrapping; first column keeps whole words; wide tables scroll.
DOC_CSS = (
    '<style>'
    'article.doc{overflow-wrap:normal}'
    'article.doc .tablewrap{margin:12px 0}'
    'article.doc table{width:max-content;max-width:100%;table-layout:auto}'
    'article.doc th,article.doc td{word-break:normal;overflow-wrap:break-word}'
    'article.doc th:first-child,article.doc td:first-child{width:1%;white-space:nowrap;overflow-wrap:normal}'
    'article.doc code{overflow-wrap:anywhere}'
    '</style>'
)


def normalize(rel: str) -> str | None:
    rel = unquote(rel or '').replace('\\', '/').lstrip('/')
    if not rel or '\0' in rel:
        return None
    parts = rel.split('/')
    if any(p in ('', '.', '..') for p in parts):
        return None
    return '/'.join(parts)


def allowed(rel: str) -> bool:
    if rel in ROOT_DOCS:
        return True
    if rel.startswith('docs/') and rel.endswith('.md'):
        return True
    if rel.startswith('docs/images/') and not rel.endswith('/'):
        return True
    return False


def resolve(rel: str):
    """Return (path, repo-relative posix) for an allowlisted file, or None."""
    rel = normalize(rel)
    if rel is None or not allowed(rel):
        return None
    root = ROOT.resolve()
    target = (ROOT / rel).resolve()
    try:
        got = target.relative_to(root).as_posix()
    except ValueError:
        return None
    if not allowed(got) or not target.is_file():
        return None
    return target, got


def _join_rel(base_dir: str, link_path: str) -> str | None:
    parts = ([p for p in base_dir.split('/') if p] if base_dir not in ('', '.') else []) + link_path.split('/')
    stack: list[str] = []
    for p in parts:
        if p in ('', '.'):
            continue
        if p == '..':
            if not stack:
                return None
            stack.pop()
            continue
        stack.append(p)
    return '/'.join(stack)


def rewrite_urls(fragment: str, doc_rel: str) -> str:
    base_dir = Path(doc_rel).parent.as_posix()

    def fix(url: str) -> str:
        url = url.strip()
        if not url or url.startswith(('#', 'mailto:', 'http://', 'https://', 'data:', '//')):
            return url
        parsed = urlparse(url)
        path = unquote(parsed.path or '')
        if path.startswith('/'):
            rel = normalize(path)
            return ('/' + rel) if rel and allowed(rel) else url
        rel = _join_rel(base_dir, path)
        if rel is None or not allowed(rel):
            return url
        out = '/' + rel
        if parsed.query:
            out += '?' + parsed.query
        if parsed.fragment:
            out += '#' + parsed.fragment
        return out

    return _ATTR.sub(lambda m: f'{m.group(1)}={m.group(2)}{fix(m.group(3))}{m.group(2)}', fragment)


def _markdown(text: str) -> str:
    import markdown
    html_body = markdown.markdown(text, extensions=['tables', 'fenced_code', 'sane_lists'], output_format='html')
    # Horizontal scroll when the label column needs its full word width.
    return re.sub(r'<table>.*?</table>', lambda m: '<div class="tablewrap">' + m.group(0) + '</div>', html_body, flags=re.S)


def render_html(target: Path, doc_rel: str) -> str:
    raw = target.read_text(encoding='utf-8')
    if target.suffix.lower() not in ('.md', '.markdown'):
        title = target.name
        body = '<pre>' + html.escape(raw) + '</pre>'
    else:
        title = target.stem.replace('_', ' ')
        m = _H1.search(raw)
        if m:
            title = m.group(1).strip()
            raw = _H1.sub('', raw, count=1).lstrip('\n')
        body = rewrite_urls(_markdown(raw), doc_rel)
    return ui.page(ROOT / doc_rel, title, DOC_CSS + f'<article class="panel doc">{body}</article>')


def handle(url_path: str):
    """If *url_path* is in the doc allowlist, return (status, body_bytes, content_type).

    Returns None when the path is outside the allowlist (caller continues with
    normal static serving). Allowlisted but missing/unsafe paths yield 404.
    """
    rel = normalize(url_path)
    if rel is None or not allowed(rel):
        return None
    hit = resolve(rel)
    if hit is None:
        return 404, None, None
    target, got = hit
    if got.startswith('docs/images/'):
        ctype = mimetypes.guess_type(target.name)[0] or 'application/octet-stream'
        if target.suffix == '.svg':
            ctype = 'image/svg+xml'
        body = target.read_bytes()
        if ctype.startswith(('text/', 'image/svg')):
            ctype += '; charset=utf-8'
        return 200, body, ctype
    return 200, render_html(target, got).encode('utf-8'), 'text/html; charset=utf-8'
