"""Optional data packages: fetch, verify and install the parts of the structure database that are published as separate downloads.

The repository contains the open-core Atlas (`layout_atlas/data/open_core.sqlite`).  The developer's full database has more: the MicroMap and BRENDA extensions, the
reference map sources and a `maps` table that lists them.  Those parts have their own licences, so they are released as separate packages (each with a notice and its
licence text inside) instead of being mixed into the code repository.  This command downloads a package, checks its SHA-256 against `data/packages.json` and extracts it
into `layout_atlas/data/`.  Nothing is overwritten silently: an existing file is moved to `<name>.bak` first, and symlinks are never touched.

  python3 project.py data list                      # packages, size, licence, installed or not
  python3 project.py data fetch micromap            # download + verify + install one package
  python3 project.py data fetch all                 # all packages
  python3 project.py data fetch all --from DIR      # install from a local directory (or a base URL / file:// URL) instead of the default release
"""
import hashlib
import json
import os
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

BASE = Path(__file__).resolve().parent
DATA = BASE / 'data'
MANIFEST = DATA / 'packages.json'


def load():
    return json.loads(MANIFEST.read_text())


def installed(pkg):
    return (DATA / pkg['marker']).exists()


def cmd_list(_args):
    m = load()
    print(f"Open Meta Maps data packages ({m['version']}) - base URL {m['base_url']}\n")
    for p in m['packages']:
        print(f"  {p['name']:15s} {p['bytes'] / 1e6:7.1f} MB  {'installed' if installed(p) else 'not installed'}")
        print(f"      {p['title']}")
        print(f"      licence / terms: {p['licence']}\n")


def _open(url_or_path, base_hint):
    s = str(url_or_path)
    if '://' in s:
        return urllib.request.urlopen(s)
    return open(s, 'rb')


def _download(pkg, source, dest):
    if source is None:
        loc = load()['base_url'] + pkg['file']
    elif '://' in source:
        loc = source.rstrip('/') + '/' + pkg['file']
    else:
        loc = str(Path(source) / pkg['file'])
    print(f"  downloading {loc}")
    h = hashlib.sha256()
    with _open(loc, None) as r, open(dest, 'wb') as w:
        total = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            w.write(chunk)
            h.update(chunk)
            total += len(chunk)
    if h.hexdigest() != pkg['sha256']:
        raise SystemExit(f"checksum mismatch for {pkg['file']}: expected {pkg['sha256']}, got {h.hexdigest()}; nothing was installed")
    print(f"  checksum ok ({total / 1e6:.1f} MB)")


def _safe_members(tf):
    root = DATA.resolve()
    for m in tf.getmembers():
        target = (DATA / m.name).resolve()
        if not str(target).startswith(str(root) + os.sep) and target != root:
            raise SystemExit(f'unsafe path in package: {m.name}')
        if m.issym() or m.islnk() or m.isdev():
            raise SystemExit(f'unexpected special file in package: {m.name}')
        yield m


def _install(pkg, archive):
    with tarfile.open(archive, 'r:gz') as tf:
        members = list(_safe_members(tf))
        backed_up = 0
        for m in members:
            if not m.isfile():
                continue
            target = DATA / m.name
            if target.is_symlink():
                raise SystemExit(f'{target} is a symlink (a developer checkout that already links to this data); not touching it')
            if target.exists():
                bak = target.with_name(target.name + '.bak')
                if bak.exists():
                    raise SystemExit(f'{bak} already exists; move it away first')
                target.rename(bak)
                backed_up += 1
        tf.extractall(DATA, members=members)
    print(f"  installed {len([m for m in members if m.isfile()])} files into {DATA.relative_to(BASE.parent)}" + (f" ({backed_up} existing file(s) kept as .bak)" if backed_up else ''))
    print(f"  notice: {DATA / pkg['notice']}")


def cmd_fetch(args):
    m = load()
    names = [a for a in args if not a.startswith('--')]
    source = None
    if '--from' in args:
        i = args.index('--from')
        source = args[i + 1]
        names = [a for a in names if a != source]
    if not names:
        raise SystemExit('usage: python3 project.py data fetch <name|all> [--from DIR_OR_URL]')
    pk = {p['name']: p for p in m['packages']}
    todo = list(pk) if names == ['all'] else names
    for n in todo:
        if n not in pk:
            raise SystemExit(f"unknown package '{n}'; known: {', '.join(pk)}")
    for n in todo:
        p = pk[n]
        print(f"[{n}] {p['title']}\n  terms: {p['licence']}")
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / p['file']
            _download(p, source, archive)
            _install(p, archive)


def main(argv):
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        return
    if argv[0] == 'list':
        cmd_list(argv[1:])
    elif argv[0] == 'fetch':
        cmd_fetch(argv[1:])
    else:
        raise SystemExit(f'unknown subcommand {argv[0]!r}; use list or fetch')


if __name__ == '__main__':
    main(sys.argv[1:])
