#!/usr/bin/env python3
"""config, species tree, intervals, SD units / loci / flanks / anchors from SEDEF (`prepare`),
lifted landings and the homolog between two lifted flanks.
"""
from __future__ import annotations

import bisect
import csv
import json
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple, Optional

class Region(NamedTuple):
    chrom: str
    start: int
    end: int

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


class Landing(NamedTuple):
    """Where a lifted interval ended up in the target (its largest landing site)."""
    chrom: str
    strand: str
    start: int
    end: int
    aligned_bp: int

    @property
    def region(self) -> Region:
        return Region(self.chrom, self.start, self.end)


class Homolog(NamedTuple):
    method: str        # "unit", "locus" or "one_flank": how it was found
    region: Region

# ------------------------------------------------------------------ settings
# How SD is grouped and anchored (prepare)
UNIT_GAP = 5_000              # SD pieces closer than this form one unit
MIN_ANCHOR_LEN = 5_000        # shorter non-SD anchors next to a unit are not used
LOCUS_GAP = 100_000           # SD pieces closer than this form one locus (fallback anchoring)
FLANK_LEN = 100_000           # length of a locus flank; also the longest unit anchor
MIN_SEQ_LEN = 1_000_000       # HAL sequences shorter than this are ignored
SAME_COPY_OVERLAP = 0.5       # two SEDEF intervals overlapping by this fraction of both
                              # lengths are the same SD copy

# How lifted flanks become a homolog
MIN_FLANK_LIFTED = 0.5        # fraction of a flank/anchor that must lift for it to count
LANDING_MERGE_GAP = 50_000    # lifted blocks closer than this are one landing site
MAX_HOMOLOG_LEN = 5_000_000   # a homolog longer than this and than
MAX_HOMOLOG_LEN_RATIO = 20    # this many times the unit length is rejected
MAX_FLANK_OVERLAP = 20_000    # lifted flanks may overlap by at most this much

# -------------------------------------------------------------------- config


class Config:
    def __init__(self, path: str) -> None:
        self.path = Path(path).resolve()
        with open(path) as fh:
            c = json.load(fh)
        base = self.path.parent

        def p(x):
            x = Path(x).expanduser()
            return x if x.is_absolute() else (base / x).resolve()

        for k in ("hal", "outdir", "tree", "species"):
            if k not in c:
                sys.exit(f"config: missing '{k}'")
        self.hal = p(c["hal"])
        self.out = p(c["outdir"])
        pre = c.get("hal_command_prefix", "")
        self.prefix = shlex.split(pre) if isinstance(pre, str) else list(pre)
        self.cache = str(c.get("hdf5_cache_bytes", 268435456))
        self.tree = Tree(parse_newick(c["tree"]))
        self.species = {}
        for code, s in c["species"].items():
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", code):
                sys.exit(f"config: species code '{code}' may only use letters, digits, _ . -")
            if code not in self.tree.parent:
                sys.exit(f"config: species '{code}' is not a leaf of the tree")
            self.species[code] = dict(genome=s["genome"], sedef=p(s["sedef"]),
                                      chrom_map=p(s["chrom_map"]) if s.get("chrom_map") else None,
                                      exclude=set(s.get("exclude_chroms", [])))
        self.codes = list(self.species)
        if len(self.codes) < 2:
            sys.exit("config: need at least two species")
        self.tag = "_".join(self.codes)
        for d in ("prep", "lifts", "logs"):
            (self.out / d).mkdir(parents=True, exist_ok=True)

    def hal_cmd(self, tool, *args):
        return self.prefix + [tool] + [str(a) for a in args]

    def prep(self, code, what):
        return self.out / "prep" / f"{code}.{what}"

    def lift(self, kind, q, t):
        suffix = {"flanks": "flank", "anchors": "anchor", "copies": "copy", "copies_all": "copy_all"}[kind]
        return self.out / "lifts" / f"{q}__TO__{t}.{suffix}.tsv"


def parse_newick(s: str):
    """'(GGO,(HSA,PTR)Hominini)Homininae;' -> nested (name, child, ...) tuples.
    Branch lengths and support values are ignored; unnamed internal nodes are
    named after their leaves (e.g. 'HSA+PTR')."""
    s = re.sub(r":[0-9.eE+-]+", "", s.strip().rstrip(";"))
    pos = 0

    def node():
        nonlocal pos
        kids = []
        if s[pos] == "(":
            pos += 1
            kids.append(node())
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            if s[pos] != ")":
                raise ValueError(f"tree: expected ')' at {pos}")
            pos += 1
        m = re.match(r"[^(),;]*", s[pos:])
        name = m.group(0).strip()
        pos += len(m.group(0))
        if not kids:
            if not name:
                raise ValueError("tree: unnamed leaf")
            return name
        if not name or re.fullmatch(r"[0-9.]+", name):
            name = "+".join(sorted(leaves(kids)))
        return (name, *kids)

    def leaves(kids):
        out = []
        for k in kids:
            out += [k] if isinstance(k, str) else leaves(k[1:])
        return out

    t = node()
    if isinstance(t, str):
        raise ValueError("tree: need at least two leaves")
    return t


class Tree:
    def __init__(self, node) -> None:
        self.parent: dict[str, str | None] = {}
        self.leaves_of: dict[str, frozenset[str]] = {}
        self.depth: dict[str, int] = {}
        self._walk(node, None, 0)

    def _walk(self, node, parent, depth) -> frozenset[str]:
        if isinstance(node, str):
            name, leaves = node, frozenset([node])
        else:
            name = node[0]
            leaves = frozenset().union(*(self._walk(c, name, depth + 1) for c in node[1:]))
        if name in self.parent:
            raise ValueError(f"tree: duplicate node name {name}")
        self.parent[name] = parent
        self.leaves_of[name] = leaves
        self.depth[name] = depth
        return leaves

    def mrca(self, present) -> str:
        present = frozenset(present)
        node = next(iter(present))
        while node is not None:
            if present <= self.leaves_of[node]:
                return node
            node = self.parent[node]
        raise RuntimeError(f"no MRCA for {sorted(present)}")


def merge_iv(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = []
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            out.append((cs, ce))
            cs, ce = s, e
    out.append((cs, ce))
    return out


def chain_iv(merged, gap):
    """(start, end, sd_bp, n_pieces) for pieces chained across gaps <= gap."""
    if not merged:
        return []
    out = []
    cs, ce = merged[0]
    sd, n = ce - cs, 1
    for s, e in merged[1:]:
        if s - ce <= gap:
            sd += e - s
            ce = e
            n += 1
        else:
            out.append((cs, ce, sd, n))
            cs, ce, sd, n = s, e, e - s, 1
    out.append((cs, ce, sd, n))
    return out


class IntervalIndex:
    """Sorted, merged intervals per chromosome with bp-overlap queries."""

    def __init__(self, by_chrom) -> None:
        self.ivs = {c: merge_iv(v) for c, v in by_chrom.items()}
        self.starts = {c: [s for s, _ in v] for c, v in self.ivs.items()}

    def overlap_bp(self, chrom, start, end) -> int:
        ivs = self.ivs.get(chrom)
        if not ivs or end <= start:
            return 0
        i = max(bisect.bisect_right(self.starts[chrom], start) - 1, 0)
        total = 0
        while i < len(ivs) and ivs[i][0] < end:
            s, e = ivs[i]
            total += max(0, min(e, end) - max(s, start))
            i += 1
        return total


def read_tsv(path: Path, ints=()):
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for k in ints:
            r[k] = int(r[k])
    return rows


def chrom_sizes(cfg: Config, code: str) -> dict[str, int]:
    path = cfg.prep(code, "sizes.tsv")
    if not path.is_file():
        cmd = cfg.hal_cmd("halStats", "--chromSizes", cfg.species[code]["genome"], cfg.hal)
        try:
            txt = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
        except (OSError, subprocess.CalledProcessError) as e:
            sys.exit(f"halStats failed for {code}: {' '.join(cmd)}\n{getattr(e, 'stderr', e)}")
        path.write_text(txt)
    sizes = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            sizes[f[0]] = int(f[1])
    return sizes


def sedef_pairs(cfg: Config, code: str):
    """SEDEF pairs in HAL coordinates -> (pairs, singles, dropped, sizes).

    pairs:   [((chrom, s, e), (chrom, s, e))] with both sides kept (families).
    singles: sides kept whose partner was dropped; still SD sequence, so they
             count towards the SD annotation but link no family.

    Columns 1-3 and 10-12 of the SEDEF bed are the two sides of each pair.
    Chromosome names are renamed through chrom_map (two columns: bed name,
    HAL name) when given; pairs on sequences not in the HAL, shorter than
    MIN_SEQ_LEN, or listed in exclude_chroms (HAL names) are dropped."""
    sizes = chrom_sizes(cfg, code)
    keep = {c for c, n in sizes.items() if n >= MIN_SEQ_LEN} - cfg.species[code]["exclude"]
    rename = None
    if cfg.species[code]["chrom_map"]:
        rename = {}
        for line in cfg.species[code]["chrom_map"].read_text().splitlines():
            f = line.split()
            if len(f) >= 2 and not line.startswith("#"):
                rename[f[0]] = f[1]
    out, singles, dropped = [], [], 0
    with cfg.species[code]["sedef"].open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                sys.exit(f"{cfg.species[code]['sedef']}: expected SEDEF bedpe-style rows "
                         "(sides in columns 1-3 and 10-12)")
            side = []
            for c, s, e in ((f[0], f[1], f[2]), (f[9], f[10], f[11])):
                c = rename.get(c) if rename is not None else c
                side.append((c, int(s), int(e)) if c in keep else None)
            if side[0] and side[1]:
                out.append(tuple(side))
            else:
                dropped += 1
                singles += [x for x in side if x]
    if not out:
        sys.exit(f"{code}: no SEDEF pair maps onto HAL genome {cfg.species[code]['genome']} "
                 f"(chromosome names differ? give a chrom_map)")
    return out, singles, dropped, sizes


def cmd_prepare(cfg: Config, args=None) -> None:
    for code in cfg.codes:
        pairs, singles, dropped, sizes = sedef_pairs(cfg, code)
        by = defaultdict(list)
        for c, s, e in [x for p in pairs for x in p] + singles:
            by[c].append((s, e))
        n_loci = n_full = n_units = 0
        with cfg.prep(code, "sd.bed").open("w") as sd_fh, \
             cfg.prep(code, "loci.tsv").open("w") as lf, \
             cfg.prep(code, "flanks.bed").open("w") as ff, \
             cfg.prep(code, "units.tsv").open("w") as ut, \
             cfg.prep(code, "anchors.bed").open("w") as ab:
            lf.write("locus_id\tchrom\tstart\tend\tsd_bp\tn_pieces\tfull_flanks\n")
            ut.write("unit_id\tchrom\tstart\tend\tsd_bp\tn_pieces\tlocus_id\tl_len\tr_len\n")
            for chrom in sorted(by):
                merged = merge_iv(by[chrom])
                size = sizes[chrom]
                for s, e in merged:
                    sd_fh.write(f"{chrom}\t{s}\t{e}\n")
                loci = chain_iv(merged, LOCUS_GAP)
                for s, e, bp, n in loci:
                    full = int(s - FLANK_LEN >= 0 and e + FLANK_LEN <= size)
                    lid = f"{code}|{chrom}|{s}|{e}"
                    lf.write(f"{lid}\t{chrom}\t{s}\t{e}\t{bp}\t{n}\t{full}\n")
                    n_loci += 1
                    if full:
                        n_full += 1
                        ff.write(f"{chrom}\t{s - FLANK_LEN}\t{s}\t{lid}|L\t0\t+\n")
                        ff.write(f"{chrom}\t{e}\t{e + FLANK_LEN}\t{lid}|R\t0\t+\n")
                units = chain_iv(merged, UNIT_GAP)
                li = 0
                for i, (s, e, bp, npc) in enumerate(units):
                    while loci[li][1] < e:
                        li += 1
                    lid = f"{code}|{chrom}|{loci[li][0]}|{loci[li][1]}"
                    prev_e = units[i - 1][1] if i else 0
                    next_s = units[i + 1][0] if i + 1 < len(units) else size
                    l0, r1 = max(s - FLANK_LEN, prev_e), min(e + FLANK_LEN, next_s)
                    l_len = s - l0 if s - l0 >= MIN_ANCHOR_LEN else 0
                    r_len = r1 - e if r1 - e >= MIN_ANCHOR_LEN else 0
                    uid = f"{code}|{chrom}|{s}|{e}"
                    ut.write(f"{uid}\t{chrom}\t{s}\t{e}\t{bp}\t{npc}\t{lid}\t{l_len}\t{r_len}\n")
                    if l_len:
                        ab.write(f"{chrom}\t{l0}\t{s}\t{uid}|L\t0\t+\n")
                    if r_len:
                        ab.write(f"{chrom}\t{e}\t{r1}\t{uid}|R\t0\t+\n")
                    n_units += 1
        sd_bp = sum(e - s for v in by.values() for s, e in merge_iv(v))
        print(f"[prepare] {code}: {len(pairs)} SEDEF pairs ({dropped} with a side on an excluded, "
              f"missing or < {MIN_SEQ_LEN // 1_000_000} Mbp sequence), {sd_bp / 1e6:.1f} Mbp SD, "
              f"{n_units} units, {n_loci} loci ({n_full} with full flanks)", flush=True)


def largest_landing(rows):
    """Largest same-chrom, same-strand cluster of lifted blocks (gaps <= LANDING_MERGE_GAP).
    rows: (chrom, start, end, strand) -> Landing."""
    best = cur = None
    for chrom, s, e, strand in sorted(rows, key=lambda r: (r[0], r[3], r[1])):
        if cur and cur[0] == chrom and cur[1] == strand and s - cur[3] <= LANDING_MERGE_GAP:
            cur[3] = max(cur[3], e)
            cur[4] += e - s
        else:
            if cur and (best is None or cur[4] > best[4]):
                best = cur
            cur = [chrom, strand, s, e, e - s]
    if cur and (best is None or cur[4] > best[4]):
        best = cur
    return Landing(*best) if best else None


def best_landing_per_name(fi, fo) -> int:
    """--noDupes halLiftover bed -> one main landing per name."""
    blocks = defaultdict(list)
    for line in fi:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6:
            blocks[f[3]].append((f[0], int(f[1]), int(f[2]), f[5]))
    for name, rows in blocks.items():
        L = largest_landing(rows)
        fo.write(f"{name}\t{L.chrom}\t{L.strand}\t{L.start}\t{L.end}\t{L.aligned_bp}\n")
    return len(blocks)


def lift_source(cfg: Config, kind: str, q: str) -> Path:
    return cfg.prep(q, {"flanks": "flanks.bed", "anchors": "anchors.bed",
                        "copies": "copies.bed", "copies_all": "copies.bed"}[kind])


def load_lifts(path: Path) -> dict[str, Landing]:
    """name -> its best landing (from a --noDupes lift)."""
    out = {}
    with path.open() as fh:
        for line in fh:
            name, chrom, strand, s, e, bp = line.rstrip("\n").split("\t")
            out[name] = Landing(chrom, strand, int(s), int(e), int(bp))
    return out


def homolog_between_flanks(left: Optional[Landing], right: Optional[Landing], span: int,
                           need_left=FLANK_LEN, need_right=FLANK_LEN) -> Optional[Region]:
    """The target interval between a lifted left and right flank, or None when
    the flanks did not lift well enough or do not land together.
    need_left/right: flank lengths; at least MIN_FLANK_LIFTED of each must lift."""
    if left is None or right is None or left.chrom != right.chrom:
        return None
    if left.aligned_bp < MIN_FLANK_LIFTED * need_left or right.aligned_bp < MIN_FLANK_LIFTED * need_right:
        return None
    if left.strand != right.strand:
        # opposite orientation: an inversion breakpoint inside the locus
        a, b = sorted((left, right), key=lambda f: f.start)
        hs, he = a.end, b.start
    elif left.strand == "+":
        hs, he = left.end, right.start
    else:
        hs, he = right.end, left.start
    if he < hs:
        if hs - he > MAX_FLANK_OVERLAP:
            return None
        hs = he = (hs + he) // 2
    if he - hs > max(MAX_HOMOLOG_LEN, MAX_HOMOLOG_LEN_RATIO * span):
        return None
    return Region(left.chrom, hs, he)


def one_flank_homolog(flank: Optional[Landing], side: str, span: int) -> Optional[Region]:
    """A span-long window next to one lifted flank (the other did not lift), or None."""
    if flank is None or flank.aligned_bp < MIN_FLANK_LIFTED * FLANK_LEN:
        return None
    if (side == "L") == (flank.strand == "+"):   # the locus continues to the right
        return Region(flank.chrom, flank.end, flank.end + span)
    return Region(flank.chrom, max(0, flank.start - span), flank.start)


def unit_homologs(cfg: Config, q: str, t: str, units, loci) -> dict[str, Homolog]:
    """unit_id -> its homolog in t: from the unit's own anchors, else its
    locus' two flanks, else one locus flank."""
    anchor = load_lifts(cfg.lift("anchors", q, t))
    flank = load_lifts(cfg.lift("flanks", q, t))
    by_locus = {}
    for l in loci:
        lid, span = l["locus_id"], l["end"] - l["start"]
        L, R = flank.get(lid + "|L"), flank.get(lid + "|R")
        h = homolog_between_flanks(L, R, span) if l["full_flanks"] else None
        if h:
            by_locus[lid] = Homolog("locus", h)
            continue
        hL, hR = one_flank_homolog(L, "L", span), one_flank_homolog(R, "R", span)
        if (hL or hR) and not (hL and hR and hL.chrom != hR.chrom):
            by_locus[lid] = Homolog("one_flank", hL or hR)
    out = {}
    for u in units:
        h = None
        if u["l_len"] and u["r_len"]:
            h = homolog_between_flanks(anchor.get(u["unit_id"] + "|L"), anchor.get(u["unit_id"] + "|R"),
                                       u["end"] - u["start"], u["l_len"], u["r_len"])
        if h:
            out[u["unit_id"]] = Homolog("unit", h)
        elif u["locus_id"] in by_locus:
            out[u["unit_id"]] = by_locus[u["locus_id"]]
    return out