#!/usr/bin/env python3
"""SD units, loci, flanks and anchors from each genome's SEDEF calls (prep/<CODE>.*).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def main():
    a = sys.argv[1:]
    cfg = C.load(a[1] if a[:1] == ["-c"] else None)
    C.O.cmd_prepare(cfg)


if __name__ == "__main__":
    main()