#!/usr/bin/env python3
"""Cactus lift of one genome's locus flanks or unit anchors to another genome."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
from common import O  # noqa: E402


def main():
    a = sys.argv[1:]
    cfgp = None
    if a and a[0] == "-c":
        cfgp, a = a[1], a[2:]
    cfg = C.load(cfgp)
    kind, q, t = a
    out = cfg.lift(kind, q, t)
    if out.is_file():
        return
    raw = cfg.out / "lifts" / f"{q}__TO__{t}.{kind}.raw.bed"
    if not C.hal_liftover(cfg, q, O.lift_source(cfg, kind, q), t, raw):
        sys.exit(f"halLiftover {kind} {q} -> {t} failed")
    with raw.open() as fi, out.with_suffix(".tmp").open("w") as fo:
        O.best_landing_per_name(fi, fo)
    out.with_suffix(".tmp").rename(out)
    raw.unlink()


if __name__ == "__main__":
    main()