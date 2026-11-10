"""Router binaries for machines that are not the arm64 Mac the checked-in binaries were built on.

`routing.py` (whose hash is recorded in saved experiment results) looks for `BASE/build/avoid_router`; the checked-in
binaries are arm64.  On another architecture this module builds its own copies from the same pinned libavoid sources
into build_native/<arch>/build/ and points `routing.BASE` there.  Neither routing.py nor the arm64 binaries change.
"""
import platform
import subprocess

from core import *
import routing


def use_native_binaries():
    arch = platform.machine()
    if arch in ('arm64', 'aarch64'):
        return None
    root = BASE / 'build_native' / arch
    build = root / 'build'
    build.mkdir(parents=True, exist_ok=True)
    for name in ('route.cpp', 'route_labels.cpp'):
        target = build / name
        if not target.exists() or target.read_bytes() != (BASE / 'build' / name).read_bytes():
            target.write_bytes((BASE / 'build' / name).read_bytes())
    link = build / 'source'
    if not link.exists():
        link.symlink_to(BASE / 'build/source')
    sources = sorted((BASE / 'build/source/cola/libavoid').glob('*.cpp'))
    for binary, cpp in (('avoid_router', 'route.cpp'), ('avoid_labels_router', 'route_labels.cpp')):
        if (build / binary).exists():
            continue
        cmd = ['clang++', '-std=c++11', '-O1', '-DNDEBUG', '-I' + str(BASE / 'build/source/cola'), str(build / cpp), *map(str, sources), '-o', str(build / binary)]
        result = subprocess.run(cmd, text=True, capture_output=True, timeout=900)
        dump(root / (binary + '.build.json'), dict(architecture=arch, command=cmd, exit_code=result.returncode, stderr=result.stderr[-4000:],
                                                    source_commit=read(OUT / 'routing_build.json')['commit']))
        if result.returncode:
            raise RuntimeError(result.stderr[-1500:])
    routing.BASE = root
    return root
