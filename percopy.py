#!/usr/bin/env python3
"""Per-copy dating of SD that has no positional call.
"""
import csv
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

SEG_LEN, SEG_MIN = 10_000, 2_000
MAX_HITS = 300
HOM_COV = 0.5
MERGE_GAP = 2_000
MAX_PER_GENOME, MAX_LEAVES = 60, 150
MIN_GROUP = 3
CN_STRONG = 3.0              # mean copies under X must also be >= this x the sister max
CN_RATIO = 1.5               # copy-number guard: mean copies under X >= this x max copies in the sister lineage
PAIR_MIN_FRAC = 0.3          # all-vs-all: aligned length >= this x shorter sequence, else the pair is "far"
FAR = 0.3                    # p-distance for pairs that do not align
ALN_TIMEOUT = 120            # s per segment for the all-vs-all alignment
NOT_STEPWISE = {"tile", "tile_in_range", "no_homolog", "range", "tile_too_old"}


def d(cfg, *p):
    return cfg.d("percopy", *p)


# ------------------------------------------------------------------ select
def select(cfg):
    import pysam
    import combine as CB
    for q in cfg.codes:
        fa = pysam.FastaFile(str(cfg.fasta[q]))
        by = defaultdict(list)
        for c, x, y, b, src, lo, hi in CB.combined_rows(cfg, q, use_copy=False):
            if src in NOT_STEPWISE and b != q:        # a species-specific tile call cannot get younger
                by[c].append((x, y))
        n = bp = 0
        with (d(cfg, "queries") / f"{q}.fa.tmp").open("w") as fh:
            for c, ivs in by.items():
                for s, e in C.O.merge_iv(ivs):
                    k = max(1, round((e - s) / SEG_LEN))
                    step = (e - s) / k
                    for i in range(k):
                        a, b = s + round(i * step), s + round((i + 1) * step)
                        if b - a < SEG_MIN:
                            continue
                        fh.write(f">{q}|{c}|{a}|{b}\n{fa.fetch(c, a, b)}\n")
                        n += 1
                        bp += b - a
        (d(cfg, "queries") / f"{q}.fa.tmp").rename(d(cfg, "queries") / f"{q}.fa")
        print(f"[percopy select] {q}: {n} segments, {bp / 1e6:.1f} Mbp", flush=True)


def map_to(cfg, t, threads):
    qs = " ".join(str(d(cfg, "queries") / f"{q}.fa") for q in cfg.codes)
    out = d(cfg, "paf") / f"to_{t}.paf"
    subprocess.run(f"cat {qs} | minimap2 -t {threads} -x asm20 -N {MAX_HITS} -p 0.2 --secondary=yes "
                   f"{cfg.fasta[t]} - > {out}.tmp 2> {out}.log && mv {out}.tmp {out}", shell=True, check=True)


# ------------------------------------------------------------------ trees
def homologs(cfg, q):
    """segment -> genome -> [(chrom, start, end)] loci covering >= HOM_COV of the segment"""
    hits = defaultdict(lambda: defaultdict(list))
    for t in cfg.codes:
        with (d(cfg, "paf") / f"to_{t}.paf").open() as fh:
            for line in fh:
                if not line.startswith(q + "|"):
                    continue
                f = line.split("\t", 12)
                hits[f[0]][(t, f[5])].append((int(f[7]), int(f[8]), int(f[2]), int(f[3])))
    out = {}
    for seg, per in hits.items():
        slen = int(seg.split("|")[3]) - int(seg.split("|")[2])
        res = defaultdict(list)
        for (t, chrom), hs in per.items():
            hs.sort()
            cur = None
            for ts, te, qs, qe in hs + [(None, None, None, None)]:
                if ts is not None and cur and ts <= cur[1] + MERGE_GAP:
                    cur[1] = max(cur[1], te)
                    cur[2].append((qs, qe))
                    continue
                if cur and sum(b - a for a, b in C.O.merge_iv(cur[2])) >= HOM_COV * slen:
                    res[t].append((chrom, cur[0], cur[1]))
                if ts is not None:
                    cur = [ts, te, [(qs, qe)]]
        out[seg] = res
    return out


def nj(names, dm):
    """neighbour joining -> unrooted tree as adjacency {node: {nbr: length}}"""
    adj = defaultdict(dict)
    act = list(names)
    dd = {a: dict(dm[a]) for a in names}
    k = 0
    while len(act) > 2:
        n = len(act)
        r = {a: sum(dd[a][b] for b in act if b != a) for a in act}
        best = None
        for i, a in enumerate(act):
            for b in act[i + 1:]:
                v = (n - 2) * dd[a][b] - r[a] - r[b]
                if best is None or v < best[0]:
                    best = (v, a, b)
        _, a, b = best
        k += 1
        u = f"#{k}"
        la = max(0.5 * dd[a][b] + (r[a] - r[b]) / (2 * (n - 2)), 0)
        adj[u][a] = adj[a][u] = la
        adj[u][b] = adj[b][u] = max(dd[a][b] - la, 0)
        dd[u] = {}
        for c in act:
            if c not in (a, b):
                dd[u][c] = dd[c][u] = 0.5 * (dd[a][c] + dd[b][c] - dd[a][b])
        act = [c for c in act if c not in (a, b)] + [u]
    a, b = act
    adj[a][b] = adj[b][a] = dd[a][b]
    return adj


class Recon:
    """species-tree helpers for LCA reconciliation"""

    def __init__(self, cfg):
        self.cfg = cfg
        self.path = {n: cfg.path_to_root(n) for n in cfg.tree.parent}
        self.depth = cfg.tree.depth

    def lca(self, sp):
        sp = list(sp)
        anc = set(self.path[sp[0]])
        for s in sp[1:]:
            anc &= set(self.path[s])
        return max(anc, key=lambda n: self.depth[n])


def reconcile(adj, leaf_sp, R):
    """root on every edge; keep the rooting with the fewest duplications.
    returns (parent map, species map M, dup set) of the best rooting"""
    edges = {tuple(sorted((a, b))) for a in adj for b in adj[a]}
    best = None
    for a, b in edges:
        parent = {a: "ROOT", b: "ROOT"}
        order = []
        stack = [(a, b), (b, a)]
        while stack:
            u, frm = stack.pop()
            order.append(u)
            for v in adj[u]:
                if v != frm:
                    parent[v] = u
                    stack.append((v, u))
        kids = defaultdict(list)
        for v, p in parent.items():
            kids[p].append(v)
        M, dup = {}, set()
        for u in reversed(order):
            if u in leaf_sp:
                M[u] = leaf_sp[u]
            else:
                M[u] = R.lca([M[k] for k in kids[u]])
                if any(M[k] == M[u] for k in kids[u]):
                    dup.add(u)
        root_m = R.lca([M[a], M[b]])
        root_dup = root_m in (M[a], M[b])
        score = len(dup) + root_dup
        if best is None or score < best[0]:
            M["ROOT"] = root_m
            best = (score, dict(parent), M, dup | ({"ROOT"} if root_dup else set()), kids)
    return best[1:]


def leaves_under(u, kids, leaf_set):
    out, stack = [], [u]
    while stack:
        x = stack.pop()
        if x in leaf_set:
            out.append(x)
        stack += kids.get(x, [])
    return out


G = {}


def fasta(t):
    import pysam
    fs = G.setdefault("fas", {})
    if t not in fs:
        fs[t] = pysam.FastaFile(str(G["cfg"].fasta[t]))
    return fs[t]


def _one(item):
    """one segment -> (row, accepted bed row | None); runs in a worker (state in G, set before the fork)"""
    cfg, q, R, idx, qcalls = G["cfg"], G["q"], G["R"], G["idx"], G["qcalls"]
    seg, per = item[:2]
    full = item[2] if len(item) > 2 else None      # rescue: full copy counts when `per` is a subsample
    _, chrom, s, e = seg.split("|")
    s, e = int(s), int(e)
    row = dict(segment=seg, chrom=chrom, start=s, end=e)
    counts = full or {t: len(per.get(t, [])) for t in cfg.codes}
    row["copies"] = ",".join(f"{t}:{counts[t]}" for t in cfg.codes if counts[t])
    # tile/current call of the segment (majority)
    cur = Counter()
    for x, y, lab in C.walk(qcalls, chrom, s, e):
        if lab:
            cur[lab] += y - x
    (cb, csrc), _ = cur.most_common(1)[0] if cur else (("undated", "none"), 0)
    row["current"] = f"{cb}:{csrc}"
    n_per = [len(v) for v in per.values()] or [0]
    if max(n_per) > MAX_PER_GENOME or sum(n_per) > MAX_LEAVES:
        row.update(status="skipped_large", gain="")
        return row, None
    leaves, leaf_sp, seqs, self_leaf = [], {}, {}, None
    for t in cfg.codes:
        for i, (c2, a, b) in enumerate(per.get(t, [])):
            name = f"{t}~{i}"
            leaves.append(name)
            leaf_sp[name] = t
            seqs[name] = (t, c2, a, b)
            if t == q and c2 == chrom and a < e and b > s:
                if self_leaf is None or abs(a - s) < abs(seqs[self_leaf][2] - s):
                    self_leaf = name
    if self_leaf is None or len(leaves) < 3 or len(set(leaf_sp.values())) < 2:
        row.update(status="no_tree", gain="")
        return row, None
    fa = Path(tempfile.mkstemp(suffix=".fa", dir=G["tmp"])[1])
    with fa.open("w") as fh:
        for name, (t, c2, a, b) in seqs.items():
            fh.write(f">{name}\n{fasta(t).fetch(c2, a, b)}\n")
    try:
        p = subprocess.run(f"minimap2 -t 1 -x asm20 -X -c -N {len(leaves) + 10} -p 0.5 {fa} {fa} 2>/dev/null",
                           shell=True, capture_output=True, text=True, timeout=ALN_TIMEOUT)
    except subprocess.TimeoutExpired:
        subprocess.run(f"pkill -f {fa}", shell=True)
        fa.unlink(missing_ok=True)
        row.update(status="timeout", gain="")
        return row, None
    fa.unlink(missing_ok=True)
    lens = {n: seqs[n][3] - seqs[n][2] for n in leaves}
    pm, pa = defaultdict(int), defaultdict(int)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if f[0] == f[5]:
            continue
        k = tuple(sorted((f[0], f[5])))
        pm[k] += int(f[9])
        pa[k] += int(f[10])
    dm = {a: {} for a in leaves}
    for i, a in enumerate(leaves):
        for b in leaves[i + 1:]:
            k = tuple(sorted((a, b)))
            dist = 1 - pm[k] / pa[k] if pa[k] >= PAIR_MIN_FRAC * min(lens[a], lens[b]) else FAR
            dm[a][b] = dm[b][a] = dist
    adj = nj(leaves, dm)
    parent, M, dup, kids = reconcile(adj, leaf_sp, R)
    # lowest duplication above the segment's own leaf
    u = parent[self_leaf]
    while u not in dup and u != "ROOT":
        u = parent[u]
    if u not in dup:
        row.update(status="no_duplication", gain="")
        return row, None
    X = M[u]
    while parent.get(u) in dup and M[parent[u]] == X:     # nested duplications on the same branch: one group
        u = parent[u]
    row["dup_branch"] = X
    under = leaves_under(u, kids, set(leaves))
    qcopies = [n for n in under if leaf_sp[n] == q]
    # copy number: genomes under X vs the sister lineage
    clade = cfg.tree.leaves_of[X] & set(cfg.codes)
    sister = (cfg.tree.leaves_of[cfg.tree.parent[X]] - cfg.tree.leaves_of[X]) & set(cfg.codes) \
        if cfg.tree.parent[X] else set()
    g_in, g_sis = sum(counts[t] for t in clade) / len(clade), max((counts[t] for t in sister), default=0)
    row["cn_genome_in_vs_sister"] = f"{g_in:.1f}/{g_sis}"
    # copy number inside the gene tree: under P, the lowest node above the duplication that holds sister copies,
    # copies of each genome under X vs separate copy groups of each sister genome (a clade of sister-lineage
    # copies only counts once, so an independent expansion in the sister lineage does not mask this one)
    P = parent.get(u)
    while P is not None and not any(leaf_sp[n] in sister for n in leaves_under(P, kids, set(leaves))):
        P = parent.get(P)
    pl = leaves_under(P, kids, set(leaves)) if P is not None else leaves
    mean_in = sum(1 for n in pl if leaf_sp[n] in clade) / len(clade)
    groups = Counter()

    def pure(n):
        lv = leaves_under(n, kids, set(leaves))
        if all(leaf_sp[x] in sister for x in lv):
            for t in {leaf_sp[x] for x in lv}:
                groups[t] += 1
            return
        for k in kids.get(n, []):
            pure(k)
    if P is not None:
        pure(P)
    max_sis = max(groups.values(), default=0)
    row["cn_in_vs_sister"] = f"{mean_in:.1f}/{max_sis}"
    # originals: other Q copies under the duplication with a stepwise call older than X
    orig = 0
    for n in qcopies:
        if n == self_leaf:
            continue
        _, c2, a, b = seqs[n]
        cc = Counter()
        for x, y, lab in C.walk(qcalls, c2, a, b):
            if lab and lab[1] == "stepwise" and lab[0] in idx and idx[lab[0]] > idx[X]:
                cc[lab[0]] += y - x
        if sum(cc.values()) >= 0.5 * (b - a):
            orig += 1
    row["q_copies_under_dup"] = len(qcopies)
    row["originals"] = orig
    st = guard(mean_in, max_sis, g_in, g_sis, orig, len(qcopies))
    row.update(status=st, gain=X if st == "gain" else "")
    return row, ((chrom, s, e, X) if st == "gain" else None)


def _cn_ok(n_in, n_sis):
    return n_in > n_sis and n_in >= CN_RATIO * n_sis


def guard(t_in, t_sis, g_in, g_sis, orig, n_q):
    """status of a duplication on branch X from genome-wide copy numbers (g_: mean copies of the genomes under X vs
    the sister genome with most): gain when the genomes under X average more copies than, at least 1.5x
    (CN_RATIO) and at least 3x (CN_STRONG) as many as the sister genome with most, and the duplication left at least
    MIN_GROUP copies in the genome. Identified originals (orig) and the in-tree count (t_) are recorded but not used:
    in a test on outgroup-dated SD, gains accepted because another copy was older were too young in 92% of old
    copies (old paralogs grouped by species in the gene tree pass as a recent duplication)."""
    if not _cn_ok(g_in, g_sis):
        return "rejected_copy_number"
    if n_q < MIN_GROUP:
        return "ambiguous_pair"
    if g_in < CN_STRONG * g_sis:
        return "weak_copy_number"
    return "gain"


def restatus(cfg, q):
    """re-apply guard() to an existing calls/<q>.tsv (no new trees) and rewrite it and the .bed"""
    f = d(cfg, "calls") / f"{q}.tsv"
    rows = list(csv.DictReader(f.open(), delimiter="\t"))
    cols = list(rows[0].keys()) if rows else []
    bed, n = [], Counter()
    for r in rows:
        if r["status"] in ("gain", "rejected_copy_number", "ambiguous_pair", "weak_no_original", "weak_copy_number"):
            t_in, t_sis = map(float, r["cn_in_vs_sister"].split("/"))
            g_in, g_sis = map(float, r["cn_genome_in_vs_sister"].split("/"))
            st = guard(t_in, t_sis, g_in, g_sis, int(r["originals"]), int(r["q_copies_under_dup"]))
            r["status"], r["gain"] = st, r["dup_branch"] if st == "gain" else ""
            if st == "gain":
                bed.append((r["chrom"], int(r["start"]), int(r["end"]), r["dup_branch"]))
        n[r["status"]] += 1
    with f.with_suffix(".tmp").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in cols) + "\n")
    with (d(cfg, "calls") / f"{q}.bed").open("w") as fh:
        for r in sorted(bed):
            fh.write("\t".join(map(str, r)) + "\n")
    f.with_suffix(".tmp").rename(f)
    print(f"[percopy restatus] {q}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)


def trees(cfg, q, threads):
    import multiprocessing as mp
    import combine as CB
    G.update(cfg=cfg, q=q, R=Recon(cfg), idx={n: i for i, n in enumerate(cfg.path_to_root(q))},
             # current (no-copy) calls of Q, to recognise originals (stepwise-dated older copies)
             qcalls=C.track([(c, x, y, (b, src)) for c, x, y, b, src, lo, hi in
                             CB.combined_rows(cfg, q, use_copy=False)]),
             tmp=tempfile.mkdtemp(dir=d(cfg, "tmp")))
    hom = homologs(cfg, q)
    print(f"[percopy trees] {q}: {len(hom)} segments with homologs, {threads} workers", flush=True)
    out_rows, bed = [], []
    with mp.get_context("fork").Pool(threads) as pool:
        for n, (row, b) in enumerate(pool.imap_unordered(_one, hom.items(), chunksize=4)):
            out_rows.append(row)
            if b:
                bed.append(b)
            if n % 1000 == 0:
                print(f"[percopy trees] {q}: {n}/{len(hom)}", flush=True)
    out_rows.sort(key=lambda r: (r["chrom"], r["start"]))
    cols = ["segment", "chrom", "start", "end", "current", "copies", "status", "gain", "dup_branch",
            "cn_in_vs_sister", "cn_genome_in_vs_sister", "q_copies_under_dup", "originals"]
    o = d(cfg, "calls") / f"{q}.tsv"
    with o.with_suffix(".tmp").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    with (d(cfg, "calls") / f"{q}.bed").open("w") as fh:
        for r in sorted(bed):
            fh.write("\t".join(map(str, r)) + "\n")
    o.with_suffix(".tmp").rename(o)
    st = Counter(r["status"] for r in out_rows)
    print(f"[percopy trees] {q}: " + ", ".join(f"{k} {v}" for k, v in st.most_common()), flush=True)




def copy_track(cfg, q):
    """accepted copy calls of Q as a track (None if the stage has not run)"""
    f = cfg.out / "percopy" / "calls" / f"{q}.bed"
    if not f.is_file():
        return None
    return C.track([(c, int(s), int(e), b) for c, s, e, b in (ln.split() for ln in f.open())])


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--config")
    ap.add_argument("cmd", choices=["select", "map", "trees", "restatus"])
    ap.add_argument("args", nargs="*")
    a = ap.parse_args()
    cfg = C.load(a.config)
    if a.cmd == "select":
        select(cfg)
    elif a.cmd == "map":
        map_to(cfg, a.args[0], int(a.args[1]) if len(a.args) > 1 else cfg.threads)
    elif a.cmd == "restatus":
        for q in (a.args or cfg.codes):
            restatus(cfg, q)
    else:
        trees(cfg, a.args[0], int(a.args[1]) if len(a.args) > 1 else 4)


if __name__ == "__main__":
    main()