#!/usr/bin/env python3
"""Shared configuration and helpers for sdplace."""
from __future__ import annotations

import bisect
import json
import os
import signal
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

import core as O

HANG_CHECK = 180        # s: empty output and no CPU use after this long = halLiftover hung opening the HAL
LIFT_TIMEOUT = 2400     # s: per halLiftover call
LIFT_TRIES = 4


class Cfg(O.Config):
    """core config plus the keys the other stages need (see config.example.json)."""

    def __init__(self, path: str) -> None:
        super().__init__(path)
        c = json.loads(Path(path).read_text())
        base = self.path.parent

        def p(x):
            x = Path(x).expanduser()
            return x if x.is_absolute() else (base / x).resolve()
        self.fasta = {}
        for code, s in c["species"].items():
            if "fasta" not in s:
                raise SystemExit(f"config: species '{code}' needs 'fasta' (the genome sequence that is in the HAL)")
            self.fasta[code] = p(s["fasta"])
        for key in ("minimap2_anchor_clades", "asm20_clades"):
            for n in c.get(key, []):
                if n not in self.tree.parent:
                    raise SystemExit(f"config: {key}: '{n}' is not a node of the tree")
        # nodes inside which unlifted rescue anchors may be placed with minimap2 (sister-level, <= ~10 My)
        self.mm2_anchor_clades = list(c.get("minimap2_anchor_clades", []))
        # nodes inside which SD is aligned with minimap2 asm20; other pairs use masked blastn. None = asm20 everywhere
        self.asm20_clades = c.get("asm20_clades")
        self.blast_priority = [x for x in c.get("blast_priority", []) if x in self.species] + \
            [x for x in self.codes if x not in c.get("blast_priority", [])]
        self.threads = int(c.get("threads", 16))
        # combine: undated runs of <= this many bp between two identical dated calls take that call (0 = off)
        self.gap_fill = int(c.get("gap_fill_max_bp", 0))
        # identity age cap for tile calls (combine.py); absent = tile calls are always used
        cap = c.get("tile_age_cap")
        self.age_cap = None
        if cap:
            ages = {n: float(v) for n, v in cap["node_ages_my"].items()}
            for n in self.tree.parent:
                if n not in self.species and n not in ages:
                    raise SystemExit(f"config: tile_age_cap.node_ages_my has no age for node '{n}'")
            ages.update({x: 0.0 for x in self.species})
            self.age_cap = dict(ages=ages, rate=float(cap.get("substitution_rate_per_site_per_my", 0.001)),
                                margin=float(cap.get("margin", 2.0)))
        self.identity_col = int(c.get("sedef_identity_column", 24))          # 1-based; SEDEF fracMatch
        self.slurm = c.get("slurm", {})
        self.root = next(n for n, par in self.tree.parent.items() if par is None)

    def d(self, *parts) -> Path:
        path = self.out.joinpath(*parts)
        path.mkdir(parents=True, exist_ok=True)
        return path

    def path_to_root(self, q):
        out = [q]
        while self.tree.parent[out[-1]]:
            out.append(self.tree.parent[out[-1]])
        return out

    def share_clade(self, q, t, nodes) -> bool:
        return any(q in self.tree.leaves_of[n] and t in self.tree.leaves_of[n] for n in nodes)

    def near(self, q, t) -> bool:
        return self.asm20_clades is None or self.share_clade(q, t, self.asm20_clades)

    def outgroup(self, node):
        """leaves of the nearest ancestor that has genomes outside `node`"""
        par = self.tree.parent[node]
        while par:
            out = (self.tree.leaves_of[par] - self.tree.leaves_of[node]) & set(self.codes)
            if out:
                return out
            par = self.tree.parent[par]
        return frozenset()


def load(path=None) -> Cfg:
    path = path or os.environ.get("SDPLACE_CONFIG")
    if not path:
        raise SystemExit("give the config: -c config.json or SDPLACE_CONFIG=config.json")
    return Cfg(path)


# ------------------------------------------------------------------ SD intervals
_SD = {}


def sd_by_chrom(cfg, code):
    """{chrom: merged SD intervals} of a genome (prep/<code>.sd.bed)"""
    if code not in _SD:
        by = defaultdict(list)
        with cfg.prep(code, "sd.bed").open() as fh:
            for line in fh:
                f = line.split("\t")
                by[f[0]].append((int(f[1]), int(f[2])))
        _SD[code] = {c: O.merge_iv(v) for c, v in by.items()}
    return _SD[code]


def sd_index(cfg, code):
    return O.IntervalIndex(sd_by_chrom(cfg, code))


def sd_pieces(by_chrom, chrom, s, e):
    return [(max(a, s), min(b, e)) for a, b in by_chrom.get(chrom, ()) if b > s and a < e]


def loci(cfg, code):
    return O.read_tsv(cfg.prep(code, "loci.tsv"), ("start", "end", "sd_bp", "n_pieces", "full_flanks"))


def units(cfg, code):
    return O.read_tsv(cfg.prep(code, "units.tsv"), ("start", "end", "sd_bp", "n_pieces", "l_len", "r_len"))


# ------------------------------------------------------------------ labelled interval tracks
def track(rows):
    """[(chrom, s, e, label)] -> {chrom: (starts, ends, labels)}, sorted, non-overlapping expected"""
    by = defaultdict(list)
    for c, s, e, lab in rows:
        by[c].append((s, e, lab))
    out = {}
    for c, v in by.items():
        v.sort()
        out[c] = ([x[0] for x in v], [x[1] for x in v], [x[2] for x in v])
    return out


def walk(tr, chrom, a, b):
    """yield (start, end, label | None) covering [a, b) from a track"""
    st, en, lab = tr.get(chrom, ([], [], []))
    i = max(bisect.bisect_right(st, a) - 1, 0)
    pos = a
    while pos < b:
        while i < len(st) and en[i] <= pos:
            i += 1
        if i < len(st) and st[i] <= pos:
            e2, l = min(b, en[i]), lab[i]
        else:
            e2, l = (min(b, st[i]) if i < len(st) else b), None
        yield pos, e2, l
        pos = e2


# ------------------------------------------------------------------ halLiftover
def _cpu_idle(pgid) -> bool:
    try:
        out = subprocess.run(["ps", "-o", "cputimes=", "-g", str(pgid)], capture_output=True, text=True).stdout
        return sum(int(x) for x in out.split()) < 5
    except Exception:
        return False


def hal_liftover(cfg, q, bed, t, out, no_dupes=True) -> bool:
    """halLiftover of `bed` (genome q) to genome t -> `out` (BED). Returns True on success.
    Runs in its own process group with no pipes, with a timeout, and retries when it hangs opening
    the HAL (empty output, no CPU): a killed wrapper otherwise leaves a hung halLiftover behind."""
    cmd = cfg.hal_cmd("halLiftover", *(["--noDupes"] if no_dupes else []), "--hdf5CacheBytes", cfg.cache, cfg.hal,
                      cfg.species[q]["genome"], bed, cfg.species[t]["genome"], out)
    out = Path(out)
    for _ in range(LIFT_TRIES):
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
        t0 = time.time()
        while p.poll() is None:
            time.sleep(5)
            el = time.time() - t0
            if el > LIFT_TIMEOUT or (el > HANG_CHECK and (not out.exists() or out.stat().st_size == 0)
                                     and _cpu_idle(p.pid)):
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                p.wait()
                break
        else:
            if p.returncode == 0:
                return True
    return False


def mbp(x) -> str:
    return f"{x / 1e6:.2f}"


__all__ = ["Cfg", "load", "sd_by_chrom", "sd_index", "sd_pieces", "loci", "units", "track", "walk", "hal_liftover",
           "mbp", "Counter", "defaultdict", "O"]