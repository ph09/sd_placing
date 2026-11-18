#!/usr/bin/env python3
"""Tile branch: date each 1 kbp tile of SD by where the tile is SD in the other genomes.
"""
import bisect
import csv
import random  # noqa: F401  (kept out of the calls: no sampling is used)
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
from common import O  # noqa: E402

TILE = 1000
MIN_TILE = 200             # a shorter remainder is joined to the previous tile
MIN_ALIGNED = 0.5
MIN_SD = 0.5
HOMOLOG_SD_BP = 1000
PRESENT, LACKING = "Ssb", "AaIin"
MAX_OG = 3                 # outgroup genomes tried per tile and round in the blastn stage
PAD = 20_000               # window padding around the unit homolog
BLAST_FLANK_MIN = 0.25     # unit homolog windows for blastn: >= this fraction of the flanks/anchors lifted
BLAST_MIN_COV = 0.25
COLS = "tile\tchrom\tstart\tend\tstatus"


def tiles_bed(cfg, q):
    return cfg.d("tiles", "bed") / f"{q}.bed"


def read_tiles(cfg, q):
    return [(c, int(s), int(e)) for c, s, e, *_ in (ln.split("\t") for ln in tiles_bed(cfg, q).open())]


def read_status(path):
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    return [(r["chrom"], int(r["start"]), int(r["end"])) for r in rows], [r["status"] for r in rows]


# ------------------------------------------------------------------ prepare / lift
def prepare(cfg):
    for q in cfg.codes:
        n = bp = 0
        with tiles_bed(cfg, q).open("w") as fh:
            for line in cfg.prep(q, "sd.bed").open():
                c, s, e = line.split("\t")[:3]
                s, e = int(s), int(e)
                pos = s
                while pos < e:
                    end = min(pos + TILE, e)
                    if e - end < MIN_TILE:
                        end = e
                    fh.write(f"{c}\t{pos}\t{end}\t{n}\t0\t+\n")
                    n += 1
                    bp += end - pos
                    pos = end
        print(f"[tiles prepare] {q}: {n} tiles, {bp / 1e6:.1f} Mbp", flush=True)


def lift(cfg, q, t):
    out = cfg.d("tiles", "lifts") / f"{q}__TO__{t}.tsv"
    if out.is_file():
        return
    raw = out.with_suffix(".raw.bed")
    if not C.hal_liftover(cfg, q, tiles_bed(cfg, q), t, raw):
        sys.exit(f"halLiftover {q}->{t} failed")
    sd = C.sd_index(cfg, t)
    aln, sdbp = Counter(), Counter()
    for line in raw.open():
        f = line.split("\t")
        aln[f[3]] += int(f[2]) - int(f[1])
        sdbp[f[3]] += sd.overlap_bp(f[0], int(f[1]), int(f[2]))
    tmp = out.with_suffix(".tmp")
    with tmp.open("w") as fh:
        for name in aln:
            fh.write(f"{name}\t{aln[name]}\t{sdbp[name]}\n")
    tmp.rename(out)
    raw.unlink()


# ------------------------------------------------------------------ Cactus status
def tile_units(cfg, q, tiles):
    units = C.units(cfg, q)
    ub = defaultdict(list)
    for u in units:
        ub[u["chrom"]].append(u)
    for v in ub.values():
        v.sort(key=lambda u: u["start"])
    ust = {c: [u["start"] for u in v] for c, v in ub.items()}

    def unit_of(c, pos):
        i = bisect.bisect_right(ust.get(c, []), pos) - 1
        return ub[c][i]["unit_id"] if i >= 0 and ub[c][i]["end"] > pos else None
    return units, [unit_of(c, (s + e) // 2) for c, s, e in tiles]


def cactus(cfg, q):
    tiles = read_tiles(cfg, q)
    units, tile_unit = tile_units(cfg, q, tiles)
    loci = C.loci(cfg, q)
    others = [t for t in cfg.codes if t != q]
    status = {}
    for t in others:
        lifted = {}
        for line in (cfg.d("tiles", "lifts") / f"{q}__TO__{t}.tsv").open():
            n, a, s = line.split("\t")
            lifted[int(n)] = (int(a), int(s))
        homs = O.unit_homologs(cfg, q, t, units, loci)
        sd_t = C.sd_index(cfg, t)
        hom_sd = {u: sd_t.overlap_bp(*h.region) for u, h in homs.items()}
        st = []
        for i, (c, s, e) in enumerate(tiles):
            a, sb = lifted.get(i, (0, 0))
            if a >= MIN_ALIGNED * (e - s):
                st.append("S" if sb >= MIN_SD * a else "A")
            else:
                u = tile_unit[i]
                st.append("I" if u in hom_sd and hom_sd[u] < HOMOLOG_SD_BP else "U")
        status[t] = st
        print(f"[tiles cactus] {q} -> {t}: " + " ".join(f"{k} {v}" for k, v in sorted(Counter(st).items())), flush=True)
    out = cfg.d("tiles", "cactus") / f"{q}.tsv"
    with out.with_suffix(".tmp").open("w") as fh:
        fh.write(COLS + "\n")
        for i, (c, s, e) in enumerate(tiles):
            fh.write(f"{i}\t{c}\t{s}\t{e}\t{''.join(status[t][i] for t in others)}\n")
    out.with_suffix(".tmp").rename(out)


# ------------------------------------------------------------------ minimap2 for undecided tiles
def queries(cfg):
    for q in cfg.codes:
        tiles, status = read_status(cfg.d("tiles", "cactus") / f"{q}.tsv")
        idx = [i for i, st in enumerate(status) if "U" in st]
        reg = cfg.d("tiles", "queries") / f"{q}.regions.txt"
        reg.write_text("".join(f"{tiles[i][0]}:{tiles[i][1] + 1}-{tiles[i][2]}\n" for i in idx))
        out = cfg.d("tiles", "queries") / f"{q}.fa"
        with out.with_suffix(".tmp").open("w") as fh:
            p = subprocess.Popen(["samtools", "faidx", str(cfg.fasta[q]), "-r", str(reg)], stdout=subprocess.PIPE, text=True)
            k = 0
            for line in p.stdout:
                if line.startswith(">"):
                    line = f">{q}|{idx[k]}\n"
                    k += 1
                fh.write(line)
            if p.wait():
                sys.exit(f"samtools faidx {q} failed")
        out.with_suffix(".tmp").rename(out)
        print(f"[tiles queries] {q}: {len(idx)} tiles with an undecided genome", flush=True)


def run_map(cfg, t, threads):
    out = cfg.d("tiles", "paf") / f"to_{t}.paf"
    qs = [str(cfg.d("tiles", "queries") / f"{c}.fa") for c in cfg.codes if c != t]
    with out.with_suffix(".tmp").open("w") as fh, out.with_suffix(".log").open("w") as log:
        cat = subprocess.Popen(["cat", *qs], stdout=subprocess.PIPE)
        rc = subprocess.call(["minimap2", "-t", str(threads), "-cx", "asm20", "--secondary=no", str(cfg.fasta[t]), "-"],
                             stdin=cat.stdout, stdout=fh, stderr=log)
        cat.wait()
    if rc:
        sys.exit(f"minimap2 to {t} exited {rc}")
    out.with_suffix(".tmp").rename(out)


def mm2calls(cfg):
    others = {q: [t for t in cfg.codes if t != q] for q in cfg.codes}
    tiles, status = {}, {}
    for q in cfg.codes:
        tiles[q], st = read_status(cfg.d("tiles", "cactus") / f"{q}.tsv")
        status[q] = [bytearray(x, "ascii") for x in st]
    for t in cfg.codes:
        sd_t = C.sd_index(cfg, t)
        for line in (cfg.d("tiles", "paf") / f"to_{t}.paf").open():
            f = line.split("\t")
            if "tp:A:P" not in f[12:]:
                continue
            ql, qs, qe = int(f[1]), int(f[2]), int(f[3])
            if qe - qs < MIN_ALIGNED * ql:
                continue
            q, idx = f[0].split("|")
            k = others[q].index(t)
            st = status[q][int(idx)]
            if st[k] not in (ord("U"), ord("s"), ord("a")):
                continue                                   # Cactus decided
            new = "s" if sd_t.overlap_bp(f[5], int(f[7]), int(f[8])) >= MIN_SD * (int(f[8]) - int(f[7])) else "a"
            if new == "s" or st[k] == ord("U"):            # several primary pieces: any on SD -> s
                st[k] = ord(new)
        print(f"[tiles mm2calls] read minimap2 hits to {t}", flush=True)
    for q in cfg.codes:
        out = cfg.d("tiles", "mm2") / f"{q}.tsv"
        with out.with_suffix(".tmp").open("w") as fh:
            fh.write(COLS + "\n")
            for i, ((c, s, e), st) in enumerate(zip(tiles[q], status[q])):
                fh.write(f"{i}\t{c}\t{s}\t{e}\t{st.decode()}\n")
        out.with_suffix(".tmp").rename(out)


# ------------------------------------------------------------------ gain rule
def make_placer(cfg):
    cache = {}

    def place(q, sts):
        """sts: genome -> status letter. -> (gain node | 'unresolved', youngest node, outgroup genomes)"""
        node = cfg.tree.mrca([q] + [t for t, x in sts.items() if x in PRESENT])
        og = cache.setdefault(node, cfg.outgroup(node))
        if og and not any(sts[t] in LACKING for t in og):
            return "unresolved", node, og
        return node, node, og
    return place


# ------------------------------------------------------------------ masked blastn in the anchored window
def _blast_job(args):
    tiles_fa, win_fa = args
    r = subprocess.run(["blastn", "-task", "blastn", "-query", tiles_fa, "-subject", win_fa, "-evalue", "1e-5",
                        "-outfmt", "6 qseqid sseqid qstart qend sstart send bitscore qlen", "-max_hsps", "50",
                        "-lcase_masking"], capture_output=True, text=True)
    return r.stdout


def blast_calls(cfg, q, t, todo, tiles, homs, sd_t, wd, threads):
    """todo: tile indices -> {tile: 's' | 'a' | 'u' (no usable hit) | None (no window)}"""
    by_win, out = defaultdict(list), {}
    for i in todo:
        if homs.get(i) is None:
            out[i] = None
        else:
            by_win[homs[i]].append(i)
    jobs = []
    for k, ((c, s, e), idx) in enumerate(by_win.items()):
        tq, tw = wd / f"{t}_{k}.q.fa", wd / f"{t}_{k}.w.fa"
        reg = "\n".join(f"{tiles[i][0]}:{tiles[i][1] + 1}-{tiles[i][2]}" for i in idx)
        fa = subprocess.run(["samtools", "faidx", str(cfg.fasta[q]), "-r", "/dev/stdin"], input=reg + "\n",
                            capture_output=True, text=True).stdout
        names = iter(idx)
        tq.write_text("".join(f">{next(names)}\n" if ln.startswith(">") else ln + "\n" for ln in fa.splitlines()))
        ws, we = max(0, s - PAD), e + PAD
        w = subprocess.run(["samtools", "faidx", str(cfg.fasta[t]), f"{c}:{ws + 1}-{we}"],
                           capture_output=True, text=True).stdout
        tw.write_text(f">{c}|{ws}\n" + w.split("\n", 1)[1] if w else "")
        jobs.append((str(tq), str(tw)))
    hits = defaultdict(list)
    with ProcessPoolExecutor(threads) as ex:
        for res in ex.map(_blast_job, jobs, chunksize=8):
            for line in res.splitlines():
                qid, sid, qs, qe, ss, se, bit, ql = line.split("\t")
                hits[int(qid)].append((float(bit), int(qs), int(qe), sid, int(ss), int(se), int(ql)))
    for i in [i for v in by_win.values() for i in v]:
        hs = sorted(hits.get(i, []), reverse=True)
        if not hs:
            out[i] = "u"                                   # no hit: absence cannot be told from divergence
            continue
        best = hs[0]
        c, ws = best[3].split("|")[0], int(best[3].split("|")[1])
        centre = (best[4] + best[5]) // 2
        cov, tbp, tsd = [], 0, 0
        for _, qs, qe, sid, ss, se, _ in hs:
            if sid == best[3] and abs((ss + se) // 2 - centre) < 5000:
                cov.append((qs, qe))
                a, z = sorted((ss, se))
                tbp += z - a + 1
                tsd += sd_t.overlap_bp(c, ws + a - 1, ws + z)
        covered = sum(e - s + 1 for s, e in O.merge_iv(cov))
        out[i] = "u" if covered < BLAST_MIN_COV * best[6] else ("s" if tsd >= MIN_SD * tbp else "a")
    for f in wd.iterdir():
        f.unlink()
    return out


def blast(cfg, q, threads):
    place = make_placer(cfg)
    others = [t for t in cfg.codes if t != q]
    tiles, status = read_status(cfg.d("tiles", "mm2") / f"{q}.tsv")
    sts = [dict(zip(others, st)) for st in status]
    units, tile_unit = tile_units(cfg, q, tiles)
    loci = C.loci(cfg, q)
    wd = Path(tempfile.mkdtemp(dir=cfg.d("tiles", "tmp")))
    hom_cache, calls = {}, {}

    def homs_for(t):
        if t not in hom_cache:
            O.MIN_FLANK_LIFTED = BLAST_FLANK_MIN
            hu = O.unit_homologs(cfg, q, t, units, loci)
            hom_cache[t] = ({i: tuple(hu[u].region) for i, u in enumerate(tile_unit) if u in hu}, C.sd_index(cfg, t))
        return hom_cache[t]

    for rnd in range(4):                       # the outgroup changes as blastn calls move the gain older
        need = defaultdict(list)
        for i, st in enumerate(sts):
            g, node, og = place(q, st)
            if g != "unresolved":
                continue
            for t in [t for t in cfg.blast_priority if t in og and st[t] == "U" and (t, i) not in calls][:MAX_OG]:
                need[t].append(i)
        print(f"[tiles blast] {q} round {rnd}: " + " ".join(f"{t}:{len(v)}" for t, v in need.items()), flush=True)
        if not need:
            break
        for t, idx in need.items():
            homs, sd_t = homs_for(t)
            res = blast_calls(cfg, q, t, idx, tiles, homs, sd_t, wd, threads)
            for i in idx:
                calls[(t, i)] = res[i]
                if res[i] in ("s", "a"):
                    sts[i][t] = {"s": "b", "a": "n"}[res[i]]
    wd.rmdir()
    out = cfg.d("tiles", "classified") / f"{q}.tsv"
    with out.with_suffix(".tmp").open("w") as fh:
        fh.write(COLS + "\tgain\tkind\toldest_possible\tyoungest_possible\n")
        for i, ((c, s, e), st) in enumerate(zip(tiles, sts)):
            g, node, og = place(q, st)
            lack = [t for t in og if st[t] in LACKING]
            kind = "" if (g == "unresolved" or not og) else ("source" if any(st[t] in "Aan" for t in lack) else "new")
            lo = cfg.tree.mrca([q] + [t for t in others if st[t] in PRESENT + "U"])
            fh.write(f"{i}\t{c}\t{s}\t{e}\t{''.join(st[t] for t in others)}\t{g}\t{kind}\t{lo}\t{node}\n")
    out.with_suffix(".tmp").rename(out)


def main():
    a = sys.argv[1:]
    cfgp = None
    if a and a[0] == "-c":
        cfgp, a = a[1], a[2:]
    cfg = C.load(cfgp)
    step, a = a[0], a[1:]
    if step == "prepare":
        prepare(cfg)
    elif step == "lift":
        lift(cfg, a[0], a[1])
    elif step == "cactus":
        cactus(cfg, a[0])
    elif step == "queries":
        queries(cfg)
    elif step == "map":
        run_map(cfg, a[0], int(a[1]) if len(a) > 1 else cfg.threads)
    elif step == "mm2calls":
        mm2calls(cfg)
    elif step == "blast":
        blast(cfg, a[0], int(a[1]) if len(a) > 1 else cfg.threads)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()