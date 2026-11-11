"""Read-only audit of generated HTML: visible UI text must not contain CJK characters.

Pathway/metabolite names inside data payloads are out of scope; this checks committed HTML shells.
"""
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CJK = re.compile(r"[\u3400-\u9fff]")
VOID = {"meta", "link", "br", "img", "input", "hr", "path", "circle", "line", "polyline"}


class P(HTMLParser):
    def __init__(self):
        super().__init__()
        self.stack = []
        self.items = []

    def handle_starttag(self, t, a):
        a = dict(a)
        if t in VOID:
            self.attr(t, a)
            return
        self.stack.append((t, a))
        self.attr(t, a)

    def attr(self, t, a):
        if self.skipped() or t == "svg":
            return
        for k in ("placeholder", "title", "aria-label", "alt"):
            if a.get(k):
                self.items.append(a.get(k))
        for k, v in a.items():
            if k.startswith("data-") and k.endswith("-zh"):
                self.items.append(f"ATTR {k}={v}")
            if k == "data-zh":
                self.items.append(f"ATTR data-zh={v}")

    def handle_endtag(self, t):
        while self.stack:
            if self.stack.pop()[0] == t:
                break

    def skipped(self):
        for t, a in self.stack:
            if t in ("script", "style", "pre", "code", "title") or "dn" in a.get("class", "").split():
                return True
        return False

    def handle_data(self, d):
        if not d.strip() or self.skipped():
            return
        if self.stack and self.stack[-1][0] == "script":
            return
        self.items.append(d)


def pages():
    yield ROOT / "index.html"
    r = ROOT / "layout_atlas/results"
    for p in sorted(r.rglob("index.html")):
        yield p


def audit(path):
    parser = P()
    parser.feed(path.read_text(encoding="utf-8"))
    return [t.strip()[:80] for t in parser.items if CJK.search(t)]


if __name__ == "__main__":
    bad = 0
    for path in pages():
        if not path.exists():
            continue
        hits = audit(path)
        bad += len(hits)
        print(f"{path.relative_to(ROOT)}: CJK-in-UI {len(hits)}")
        if "-v" in sys.argv:
            for x in hits[:20]:
                print("  ", x)
    print("total flags", bad)
    sys.exit(1 if bad else 0)
