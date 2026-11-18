#!/usr/bin/env python3
"""Stepwise branch: date each SD base inside flank-defined homologous loci.
"""
import bisect
import csv
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
from common import O  # noqa: E402

FLANK = O.FLANK_LEN
STEP = 100_000              # outer rescue window
N_WALK = 9                  # windows 1..9 beyond the flank: up to 1 Mbp out
MIN_MAPQ = 30
CHUNK = 400                 # anchors per halLiftover call
MIN_SD_OVERLAP = 1000       # a homolog with less SD than this is 'locus there, no SD'
MIN_HIT_BP = 300
HOM_COLS = "species\tlocus_id\ttarget\tq_start\tq_end\th_chrom\th_start\th_end\th_sd_bp\tlevel\tsource"


def S(cfg, *parts):
    return cfg.d("stepwise", *parts)


def good(v, length):
    return v is not None and v.aligned_bp >= O.MIN_FLANK_LIFTED * length


def flank_lifts(cfg, q):
    return {t: O.load_lifts(cfg.lift("flanks", q, t)) for t in cfg.codes if t != q}


def locus_homolog(lf, loc):
    lid = loc["locus_id"]
    return O.homolog_between_flanks(lf.get(lid + "|L"), lf.get(lid + "|R"), loc["end"] - loc["start"])


def units_of(cfg, q, loc):
    return O.chain_iv(C.sd_pieces(C.sd_by_chrom(cfg, q), loc["chrom"], loc["start"], loc["end"]), O.UNIT_GAP)


# ------------------------------------------------------------------ locus homologs
def homologs(cfg):
    sd = {c: C.sd_index(cfg, c) for c in cfg.codes}
    out = S(cfg) / "homologs_locus.tsv"
    with out.with_suffix(".tmp").open("w") as fh:
        fh.write(HOM_COLS + "\n")
        for q in cfg.codes:
            lifts = flank_lifts(cfg, q)
            n = 0
            for loc in C.loci(cfg, q):
                if not loc["full_flanks"]:
                    continue
                for t, lf in lifts.items():
                    h = locus_homolog(lf, loc)
                    if h is not None:
                        fh.write(f"{q}\t{loc['locus_id']}\t{t}\t{loc['start']}\t{loc['end']}\t{h.chrom}\t{h.start}\t"
                                 f"{h.end}\t{sd[t].overlap_bp(*h)}\tlocus\tcactus\n")
                        n += 1
            print(f"[homologs] {q}: {n} locus x genome homologs", flush=True)
    out.with_suffix(".tmp").rename(out)


def read_homologs(cfg, q=None):
    """(q, locus_id) -> t -> [(q_start, q_end, (chrom, start, end), h_sd_bp)]"""
    by = defaultdict(lambda: defaultdict(list))
    for name in ("homologs_locus.tsv", "homologs_rescue.tsv"):
        p = S(cfg) / name
        if not p.is_file():
            continue
        with p.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if q is None or r["species"] == q:
                    by[(r["species"], r["locus_id"])][r["target"]].append(
                        (int(r["q_start"]), int(r["q_end"]), (r["h_chrom"], int(r["h_start"]), int(r["h_end"])),
                         int(r["h_sd_bp"])))
    return by


# ------------------------------------------------------------------ rescue: anchors
def failing(cfg, q):
    """{T: [loci of Q with both flanks but no homolog in T]}"""
    lifts = flank_lifts(cfg, q)
    out = defaultdict(list)
    for loc in C.loci(cfg, q):
        if loc["full_flanks"]:
            for t, lf in lifts.items():
                if locus_homolog(lf, loc) is None:
                    out[t].append(loc)
    return out


def anchors(cfg):
    for q in cfg.codes:
        sizes = O.chrom_sizes(cfg, q)
        fail = failing(cfg, q)
        for t in cfg.codes:
            if t == q:
                continue
            with (S(cfg, "rescue", "anchors") / f"{q}__TO__{t}.bed").open("w") as fh:
                for loc in fail.get(t, []):
                    c, s, e, lid = loc["chrom"], loc["start"], loc["end"], loc["locus_id"]
                    # L0 / R0 = the locus flanks: already lifted with Cactus, listed for the minimap2 fallback
                    fh.write(f"{c}\t{s - FLANK}\t{s}\t{lid}|L0\t0\t+\n{c}\t{e}\t{e + FLANK}\t{lid}|R0\t0\t+\n")
                    for k in range(1, N_WALK + 1):
                        if s - FLANK - k * STEP >= 0:
                            fh.write(f"{c}\t{s - FLANK - k * STEP}\t{s - FLANK - (k - 1) * STEP}\t{lid}|L{k}\t0\t+\n")
                        if e + FLANK + k * STEP <= sizes[c]:
                            fh.write(f"{c}\t{e + FLANK + (k - 1) * STEP}\t{e + FLANK + k * STEP}\t{lid}|R{k}\t0\t+\n")
                    us = units_of(cfg, q, loc)
                    for i in range(1, len(us)):
                        fh.write(f"{c}\t{us[i - 1][1]}\t{us[i][0]}\t{lid}|G{i}\t0\t+\n")
        print(f"[anchors] {q}: " + " ".join(f"{t}={len(v)}" for t, v in fail.items()), flush=True)


def anchor_spans(cfg, q, t):
    out = {}
    with (S(cfg, "rescue", "anchors") / f"{q}__TO__{t}.bed").open() as fh:
        for line in fh:
            c, a, b, n = line.split("\t")[:4]
            out[n] = (c, int(a), int(b))
    return out


# ------------------------------------------------------------------ rescue: Cactus lifts
def rescue_lift(cfg, q, cpus):
    """Chunked, resumable (a finished chunk leaves <chunk>.out.done). Temporary files stay under <outdir> so a
    containerised halLiftover can see them."""
    tmp = S(cfg, "rescue", "tmp", q)
    targets = [t for t in cfg.codes if t != q and not (S(cfg, "rescue", "cactus") / f"{q}__TO__{t}.tsv").is_file()]
    chunks = defaultdict(list)
    for t in targets:
        lines = (S(cfg, "rescue", "anchors") / f"{q}__TO__{t}.bed").read_text().splitlines(keepends=True)
        lines = [x for x in lines if not x.split("\t")[3].endswith(("|L0", "|R0"))]
        for i in range(0, len(lines), CHUNK):
            f = tmp / f"{t}.c{i // CHUNK:04d}.bed"
            f.write_text("".join(lines[i:i + CHUNK]))
            chunks[t].append(f)

    def one(args):
        t, bed = args
        out = bed.with_suffix(".out")
        done = Path(str(out) + ".done")
        if done.is_file():
            return None
        if C.hal_liftover(cfg, q, bed, t, out):
            done.touch()
            return None
        return f"FAILED {bed} -> {t}"
    with ThreadPoolExecutor(cpus) as ex:
        errs = [e for e in ex.map(one, [(t, c) for t in targets for c in chunks[t]]) if e]
    if errs:
        sys.exit("\n".join(errs))
    for t in targets:
        out = S(cfg, "rescue", "cactus") / f"{q}__TO__{t}.tsv"

        def rows():
            for c in chunks[t]:
                yield from c.with_suffix(".out").open()
        with out.with_suffix(".tmp").open("w") as fo:
            O.best_landing_per_name(rows(), fo)
        out.with_suffix(".tmp").rename(out)
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    print(f"[rescue-lift] {q}: {len(targets)} targets", flush=True)


def cactus_anchor_lifts(cfg, q, t):
    """anchor name -> Landing: the locus flanks (as L0 / R0) and the rescue anchors"""
    out = {n[:-2] + "|" + n[-1] + "0": v for n, v in O.load_lifts(cfg.lift("flanks", q, t)).items()}
    p = S(cfg, "rescue", "cactus") / f"{q}__TO__{t}.tsv"
    if p.is_file():
        out.update(O.load_lifts(p))
    return out


# ------------------------------------------------------------------ rescue: minimap2 fallback
def rescue_mm2(cfg, t, threads):
    import pysam
    qs = [q for q in cfg.codes if q != t and cfg.share_clade(q, t, cfg.mm2_anchor_clades)
          and (S(cfg, "rescue", "cactus") / f"{q}__TO__{t}.tsv").is_file()
          and not (S(cfg, "rescue", "mm2") / f"{q}__TO__{t}.tsv").is_file()]
    tmp = Path(tempfile.mkdtemp(dir=S(cfg, "rescue", "tmp")))
    fa_q, paf = tmp / "q.fa", tmp / "out.paf"
    n = 0
    with fa_q.open("w") as fo:
        for q in qs:
            anc = anchor_spans(cfg, q, t)
            cl = cactus_anchor_lifts(cfg, q, t)
            fq = pysam.FastaFile(str(cfg.fasta[q]))
            first = {}            # outer windows beyond the first Cactus-lifted one on a side are never used
            for name, (c, a, b) in anc.items():
                lid, tg = name.rsplit("|", 1)
                if tg[0] in "LR" and good(cl.get(name), b - a):
                    first[(lid, tg[0])] = min(first.get((lid, tg[0]), 99), int(tg[1:]))
            for name, (c, a, b) in anc.items():
                lid, tg = name.rsplit("|", 1)
                if good(cl.get(name), b - a) or (tg[0] in "LR" and int(tg[1:]) > first.get((lid, tg[0]), 99)):
                    continue
                fo.write(f">{q}#{name}\n{fq.fetch(c, a, b)}\n")
                n += 1
    rows = defaultdict(lambda: defaultdict(list))
    if n:
        subprocess.check_call(f"minimap2 -t {threads} -x asm20 -c --secondary=no {cfg.fasta[t]} {fa_q} > {paf}",
                              shell=True, stderr=subprocess.DEVNULL)
        with paf.open() as fh:
            for line in fh:
                f = line.split("\t")
                if int(f[11]) >= MIN_MAPQ:
                    q, name = f[0].split("#", 1)
                    rows[q][name].append((f[5], int(f[7]), int(f[8]), f[4]))
        paf.unlink()
    for q in qs:
        out = S(cfg, "rescue", "mm2") / f"{q}__TO__{t}.tsv"
        with out.with_suffix(".tmp").open("w") as fo:
            for name, rr in rows[q].items():
                L = O.largest_landing(rr)
                fo.write(f"{name}\t{L.chrom}\t{L.strand}\t{L.start}\t{L.end}\t{L.aligned_bp}\n")
        out.with_suffix(".tmp").rename(out)
    fa_q.unlink()
    tmp.rmdir()
    print(f"[rescue-mm2] -> {t}: {n} anchors from {len(qs)} genomes", flush=True)


# ------------------------------------------------------------------ rescue: homologs
def rescue_homologs(cfg):
    sd = {c: C.sd_index(cfg, c) for c in cfg.codes}
    summ = defaultdict(Counter)
    out = S(cfg) / "homologs_rescue.tsv"
    with out.with_suffix(".tmp").open("w") as fh:
        fh.write(HOM_COLS + "\n")
        for q in cfg.codes:
            fail = failing(cfg, q)
            for t in cfg.codes:
                if t == q:
                    continue
                anc = anchor_spans(cfg, q, t)
                cl = cactus_anchor_lifts(cfg, q, t)
                p = S(cfg, "rescue", "mm2") / f"{q}__TO__{t}.tsv"
                mm = O.load_lifts(p) if p.is_file() else {}

                def pick(name, length=None):
                    if length is None:
                        if name not in anc:
                            return None, "", 0
                        length = anc[name][2] - anc[name][1]
                    if good(cl.get(name), length):
                        return cl[name], "cactus", length
                    if good(mm.get(name), length):
                        return mm[name], "mm2", length
                    return None, "", length

                for loc in fail.get(t, []):
                    lid, s, e = loc["locus_id"], loc["start"], loc["end"]
                    side = {}
                    for sd_ in "LR":
                        side[sd_] = (None, 0, "")
                        for k in range(N_WALK + 1):
                            v, src, _ = pick(f"{lid}|{sd_}{k}", FLANK if k == 0 else None)
                            if v is not None:
                                side[sd_] = (v, k * STEP, src)
                                break
                    (L, dl, sl), (R, dr, sr) = side["L"], side["R"]
                    h = O.homolog_between_flanks(L, R, e - s + dl + dr, STEP, STEP) if L and R else None
                    if h is not None:
                        src = "+".join(sorted({sl, sr}))
                        fh.write(f"{q}\t{lid}\t{t}\t{s}\t{e}\t{h.chrom}\t{h.start}\t{h.end}\t{sd[t].overlap_bp(*h)}"
                                 f"\trescue_locus\t{src}\n")
                        summ[q][f"locus:{src}"] += loc["sd_bp"]
                        continue
                    us = units_of(cfg, q, loc)
                    for i, (us_, ue, bp, _) in enumerate(us):
                        la, s1, ll = pick(f"{lid}|G{i}") if i else (L, sl, FLANK)
                        ra, s2, rl = pick(f"{lid}|G{i + 1}") if i + 1 < len(us) else (R, sr, FLANK)
                        if la is None or ra is None:
                            continue
                        h = O.homolog_between_flanks(la, ra, ue - us_, ll, rl)
                        if h is None:
                            continue
                        src = "+".join(sorted({s1, s2}))
                        fh.write(f"{q}\t{lid}\t{t}\t{us_}\t{ue}\t{h.chrom}\t{h.start}\t{h.end}\t{sd[t].overlap_bp(*h)}"
                                 f"\trescue_unit\t{src}\n")
                        summ[q][f"unit:{src}"] += bp
            print(f"[rescue-homologs] {q}", flush=True)
    out.with_suffix(".tmp").rename(out)
    keys = sorted({k for v in summ.values() for k in v})
    with (S(cfg) / "homologs_rescue_summary.tsv").open("w") as fh:
        fh.write("genome\t" + "\t".join(f"{k}_sd_mbp_x_targets" for k in keys) + "\n")
        for q in cfg.codes:
            fh.write(q + "\t" + "\t".join(C.mbp(summ[q][k]) for k in keys) + "\n")


# ------------------------------------------------------------------ per-base alignment
def _blast_cover(seqs, target, wd, tag):
    qf, sf = wd / f"{tag}.q.fa", wd / f"{tag}.s.fa"
    qf.write_text("".join(f">{i}\n{s}\n" for i, (_, s) in enumerate(seqs)))
    sf.write_text(f">t\n{target}\n")
    out = subprocess.run(["blastn", "-task", "dc-megablast", "-query", str(qf), "-subject", str(sf), "-lcase_masking",
                          "-evalue", "1e-10", "-max_hsps", "2000", "-outfmt", "6 qseqid qstart qend"],
                         capture_output=True, text=True).stdout
    iv = defaultdict(list)
    for line in out.splitlines():
        i, a, b = line.split("\t")
        a, b = sorted((int(a), int(b)))
        if b - a + 1 >= MIN_HIT_BP:
            iv[int(i)].append((a - 1, b))
    return iv


def _align_one(args):
    import mappy
    import pysam
    cfgp, q, t, near, jobs, tmpdir = args
    cfg = C.load(cfgp)
    fq, ft = pysam.FastaFile(str(cfg.fasta[q])), pysam.FastaFile(str(cfg.fasta[t]))
    sq, st = C.sd_by_chrom(cfg, q), C.sd_by_chrom(cfg, t)
    wd = Path(tempfile.mkdtemp(dir=tmpdir))
    out = []
    for n, (lid, (qc, qs, qe), (tc, ts, te)) in enumerate(jobs):
        qp, tp = C.sd_pieces(sq, qc, qs, qe), C.sd_pieces(st, tc, ts, te)
        if not qp or not tp:
            continue
        target = ("N" * 200).join(ft.fetch(tc, a, b) for a, b in tp)
        seqs = [((a, b), fq.fetch(qc, a, b)) for a, b in qp]
        if near:
            al = mappy.Aligner(seq=target, preset="asm20", n_threads=1)
            hits = {i: [(h.q_st, h.q_en) for h in al.map(s) if h.q_en - h.q_st >= MIN_HIT_BP]
                    for i, (_, s) in enumerate(seqs)}
        else:
            hits = _blast_cover(seqs, target, wd, f"{t}_{n}")
        for i, ivs in hits.items():
            a0 = seqs[i][0][0]
            for x, y in O.merge_iv(ivs):
                out.append((lid, t, qc, a0 + x, a0 + y))
    for f in wd.iterdir():
        f.unlink()
    wd.rmdir()
    return out


def align(cfg, q, threads):
    chrom = {l["locus_id"]: l["chrom"] for l in C.loci(cfg, q)}
    per_t = defaultdict(list)
    for (_, lid), d in read_homologs(cfg, q).items():
        for t, regs in d.items():
            for qs, qe, h, hsd in regs:
                if hsd >= MIN_SD_OVERLAP:
                    per_t[t].append((lid, (chrom[lid], qs, qe), h))
    tmpdir = str(S(cfg, "tmp"))
    tasks = []
    for t, jobs in per_t.items():
        jobs.sort(key=lambda j: -(j[1][2] - j[1][1]))
        near = cfg.near(q, t)
        k = threads if near else 4 * threads                 # blastn tasks are slower: spread finer
        tasks += [(str(cfg.path), q, t, near, jobs[i::k], tmpdir) for i in range(k) if jobs[i::k]]
    tasks.sort(key=lambda x: (x[3], -len(x[4])))            # blastn first
    out = S(cfg, "aln") / f"{q}.tsv"
    with ProcessPoolExecutor(threads) as ex, out.with_suffix(".tmp").open("w") as fh:
        fh.write("locus_id\ttarget\tchrom\tstart\tend\n")
        for res in ex.map(_align_one, tasks, chunksize=1):
            for r in res:
                fh.write("\t".join(map(str, r)) + "\n")
    out.with_suffix(".tmp").rename(out)
    print(f"[align] {q}: {sum(len(v) for v in per_t.values())} region x genome pairs", flush=True)


# ------------------------------------------------------------------ stepwise placement
def place(cfg):
    tree = cfg.tree
    prof = defaultdict(Counter)
    for q in cfg.codes:
        hom = read_homologs(cfg, q)
        aln = defaultdict(lambda: defaultdict(list))
        with (S(cfg, "aln") / f"{q}.tsv").open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                aln[r["locus_id"]][r["target"]].append((int(r["start"]), int(r["end"])))
        sd = C.sd_by_chrom(cfg, q)
        path = cfg.path_to_root(q)
        sister = [(tree.leaves_of[node] - tree.leaves_of[child]) & set(cfg.codes) for child, node in zip(path, path[1:])]
        out = S(cfg, "placed") / f"{q}.tsv"
        with out.with_suffix(".tmp").open("w") as fo:
            fo.write("chrom\tstart\tend\tgain\tyoungest\toldest\n")
            for loc in C.loci(cfg, q):
                lid, c = loc["locus_id"], loc["chrom"]
                regs = {t: [(qs, qe) for qs, qe, _, _ in v] for t, v in hom.get((q, lid), {}).items()}
                ivs = {t: O.merge_iv(v) for t, v in aln.get(lid, {}).items()}
                starts = {t: [x for x, _ in v] for t, v in ivs.items()}
                cuts = sorted({x for v in ivs.values() for iv in v for x in iv} |
                              {x for v in regs.values() for iv in v for x in iv})
                for a, b in C.sd_pieces(sd, c, loc["start"], loc["end"]):
                    pts = [a] + cuts[bisect.bisect_right(cuts, a):bisect.bisect_left(cuts, b)] + [b]
                    for x, y in zip(pts, pts[1:]):
                        if y <= x:
                            continue
                        pres = set()
                        for t, v in ivs.items():
                            i = bisect.bisect_right(starts[t], x) - 1
                            if i >= 0 and v[i][1] >= y:
                                pres.add(t)
                        has = {t for t, v in regs.items() if any(qs <= x and y <= qe for qs, qe in v)}
                        state = ["p" if s & pres else "l" if s & has else "m" for s in sister]
                        if all(s == "m" for s in state):
                            lab, lo, hi = "no_homolog", "", ""
                        else:
                            k = max((i for i, s in enumerate(state) if s == "p"), default=-1)
                            gain = path[k + 1] if k >= 0 else q
                            j = k + 1
                            while j < len(state) and state[j] == "m":
                                j += 1
                            if j == k + 1:                     # the node just above is lacking, or gain is the root
                                lab, lo, hi = gain, gain, gain
                            else:                              # only missing nodes up to path[j]
                                lab, lo, hi = "unresolved", gain, path[j]
                        fo.write(f"{c}\t{x}\t{y}\t{lab}\t{lo}\t{hi}\n")
                        prof[q][lab if lab != "unresolved" else f"unresolved>={lo}"] += y - x
        out.with_suffix(".tmp").rename(out)
        print(f"[place] {q}", flush=True)
    with (S(cfg) / "stepwise_genome.tsv").open("w") as fh:
        fh.write("genome\tsd_mbp\tspecies_specific_mbp\tresolved_mbp\tunresolved_mbp\tno_homolog_mbp\n")
        for q in cfg.codes:
            p = prof[q]
            tot = sum(p.values())
            un = sum(v for k, v in p.items() if k.startswith("unresolved"))
            fh.write(f"{q}\t{C.mbp(tot)}\t{C.mbp(p[q])}\t{C.mbp(tot - un - p['no_homolog'])}\t{C.mbp(un)}"
                     f"\t{C.mbp(p['no_homolog'])}\n")


def main():
    a = sys.argv[1:]
    cfgp = None
    if a and a[0] == "-c":
        cfgp, a = a[1], a[2:]
    cfg = C.load(cfgp)
    step, a = a[0], a[1:]
    n = int(a[1]) if len(a) > 1 else cfg.threads
    if step == "homologs":
        homologs(cfg)
    elif step == "anchors":
        anchors(cfg)
    elif step == "rescue-lift":
        rescue_lift(cfg, a[0], n)
    elif step == "rescue-mm2":
        rescue_mm2(cfg, a[0], n)
    elif step == "rescue-homologs":
        rescue_homologs(cfg)
    elif step == "align":
        align(cfg, a[0], n)
    elif step == "place":
        place(cfg)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()