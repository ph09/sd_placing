#!/usr/bin/env python3
"""Combined SD age. It uses the stepwise call where it resolves a branch, and the tile call where it cannot.
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def stepwise_track(cfg, q):
    rows = []
    with (cfg.d("stepwise", "placed") / f"{q}.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lab = f"U|{r['youngest']}|{r['oldest']}" if r["gain"] == "unresolved" else r["gain"]
            rows.append((r["chrom"], int(r["start"]), int(r["end"]), lab))
    return C.track(rows)


def tile_track(cfg, q):
    rows = []
    with (cfg.d("tiles", "classified") / f"{q}.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append((r["chrom"], int(r["start"]), int(r["end"]), r["gain"]))
    return C.track(rows)


def divergence_track(cfg, q):
    """Per SD base of Q: the largest divergence (1 - SEDEF identity, column sedef_identity_column) among the
    SEDEF pairs covering it, i.e. its most diverged paralog. {} when the config has no tile_age_cap."""
    if cfg.age_cap is None:
        return {}
    import heapq
    col = cfg.identity_col - 1
    rename = None
    if cfg.species[q].get("chrom_map"):
        rename = dict(ln.split()[:2] for ln in cfg.species[q]["chrom_map"].read_text().splitlines()
                      if ln.strip() and not ln.startswith("#"))
    by = C.defaultdict(list)
    with open(cfg.species[q]["sedef"]) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            try:
                d = 1.0 - float(f[col])
            except (IndexError, ValueError):
                continue
            for c, a, b in ((f[0], f[1], f[2]), (f[9], f[10], f[11])):
                c = rename.get(c, c) if rename else c
                if c:
                    by[c].append((int(a), int(b), d))
    rows = []
    for c, ivs in by.items():                        # sweep: max d over overlapping intervals
        ev = sorted({x for a, b, _ in ivs for x in (a, b)})
        ivs.sort()
        heap, i = [], 0
        for x, y in zip(ev, ev[1:]):
            while i < len(ivs) and ivs[i][0] <= x:
                heapq.heappush(heap, (-ivs[i][2], ivs[i][1]))
                i += 1
            while heap and heap[0][1] <= x:
                heapq.heappop(heap)
            if heap:
                rows.append((c, x, y, -heap[0][0]))
    return C.track(rows)


def too_old(cfg, branch, d):
    """True when the whole branch is older than divergence d allows: the youngest end of the branch (age of
    its node; 0 for a terminal branch) exceeds margin x d / (2 x rate) My."""
    cap = cfg.age_cap
    return cap["ages"][branch] > cap["margin"] * d / (2 * cap["rate"])


def combined_rows(cfg, q, use_copy=True):
    """_copy_rows, then gap fill (cfg.gap_fill > 0, only with use_copy): a run of undated bases of at most
    cfg.gap_fill bp that touches the same dated call on both sides takes that call, source 'gap_fill'. In a holdout
    of dated sequence the flanking call matched 98% (2 kbp windows) / 95% (10 kbp) of the bases between."""
    rows = _copy_rows(cfg, q, use_copy)
    if not (use_copy and cfg.gap_fill):
        yield from rows
        return
    buf = []

    def flush():
        i = 0
        while i < len(buf):
            if buf[i][3] != "undated":
                yield buf[i]
                i += 1
                continue
            j = i
            while j + 1 < len(buf) and buf[j + 1][3] == "undated" and buf[j + 1][1] == buf[j][2]:
                j += 1
            L, R = (buf[i - 1] if i else None), (buf[j + 1] if j + 1 < len(buf) else None)
            size = buf[j][2] - buf[i][1]
            if (L and R and L[2] == buf[i][1] and R[1] == buf[j][2] and L[3] == R[3] and size <= cfg.gap_fill):
                yield buf[i][0], buf[i][1], buf[j][2], L[3], "gap_fill", L[3], L[3]
            else:
                yield from buf[i:j + 1]
            i = j + 1
    for r in rows:
        if buf and (r[0] != buf[-1][0] or r[1] != buf[-1][2]):     # new chromosome or a gap: flush the block
            yield from flush()
            buf = []
        buf.append(r)
    yield from flush()


def _copy_rows(cfg, q, use_copy=True):
    """Base calls (_base_rows) with the per-copy calls of percopy.py in place of tile calls and undated SD where
    they are accepted: a copy call replaces a tile call only when it is younger, an undated base only within its
    stepwise range (youngest..oldest). source 'copy'."""
    ct = None
    if use_copy:
        import percopy
        ct = percopy.copy_track(cfg, q)
    idx = {n: i for i, n in enumerate(cfg.path_to_root(q))}
    for c, x, y, b, src, lo, hi in _base_rows(cfg, q):
        shown = (lo, hi) if b == "undated" else (b, b)      # dated rows are written with youngest = oldest = branch
        if src == "stepwise" or ct is None:
            yield c, x, y, b, src, *shown
            continue
        for x2, y2, X in C.walk(ct, c, x, y):
            ok = X is not None and X in idx and idx[lo] <= idx[X] <= idx[hi] and (b == "undated" or idx[X] < idx[b])
            if ok:
                yield c, x2, y2, X, "copy", X, X
            else:
                yield c, x2, y2, b, src, *shown


def _base_rows(cfg, q):
    """yield (chrom, start, end, branch | 'undated', source, youngest, oldest) over Q's SD.
    youngest/oldest: the range the stepwise walk allows (for a tile call too: _copy_rows needs it to place a
    younger per-copy call).
    Undated sources: no_homolog, range, tile_too_old (a tile call older than the base's most diverged SEDEF paralog
    allows, see too_old: Cactus can align a recent copy to an older paralog, which dates it too old)."""
    path = cfg.path_to_root(q)
    idx = {n: i for i, n in enumerate(path)}
    st, tt, dv = stepwise_track(cfg, q), tile_track(cfg, q), divergence_track(cfg, q)
    for c, ivs in C.sd_by_chrom(cfg, q).items():
        for a, b in ivs:
            for x, y, s in C.walk(st, c, a, b):
                if s in idx:
                    yield c, x, y, s, "stepwise", s, s
                    continue
                if s is None or s == "no_homolog":
                    lo, hi, src, und = q, path[-1], "tile", "no_homolog"
                else:
                    _, lo, hi = s.split("|")
                    src, und = "tile_in_range", "range"
                for x2, y2, t in C.walk(tt, c, x, y):
                    if t in idx and idx[lo] <= idx[t] <= idx[hi]:
                        for x3, y3, d in C.walk(dv, c, x2, y2):
                            if d is not None and too_old(cfg, t, d):
                                yield c, x3, y3, "undated", "tile_too_old", lo, hi
                            else:
                                yield c, x3, y3, t, src, lo, hi
                    else:
                        yield c, x2, y2, "undated", und, lo, hi


def flag_track(cfg, q):
    """species.<q>.flag_regions BED (merged) as a track, or None. Calls there are kept and reported separately."""
    f = cfg.flag_regions.get(q)
    if not f:
        return None
    by = defaultdict(list)
    with open(f) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            x = line.split("\t")
            by[x[0]].append((int(x[1]), int(x[2])))
    return C.track([(c, a, b, True) for c, v in by.items() for a, b in C.O.merge_iv(v)])


def flagged_bp(ft, c, a, b):
    return sum(y - x for x, y, l in C.walk(ft, c, a, b) if l) if ft else 0


UNIT_GAPS = {"units_5kb": 5_000, "units_100kb": 100_000}
UNIT_MIN_BP = 1_000


def unit_counts(cfg, q):
    """per unit kind: Counter(branch -> units), Counter(branch -> SD bp of those units)"""
    rows = C.track([(c, x, y, b) for c, x, y, b in
                    ((r["chrom"], int(r["start"]), int(r["end"]), r["branch"]) for r in
                     csv.DictReader((cfg.d("results", "combined") / f"{q}.tsv").open(), delimiter="\t"))
                    if b != "undated"])
    n, size = defaultdict(Counter), defaultdict(Counter)
    for c, ivs in C.sd_by_chrom(cfg, q).items():
        for kind, gap in UNIT_GAPS.items():
            for us, ue, sdbp, _ in C.O.chain_iv(ivs, gap):
                per = Counter()
                for x, y, b in C.walk(rows, c, us, ue):
                    if b is not None:
                        per[b] += y - x
                if per:
                    b, bp = per.most_common(1)[0]
                    if bp >= UNIT_MIN_BP:
                        n[kind][b] += 1
                        size[kind][b] += sdbp
    return n, size


def combine(cfg):
    tree = cfg.tree
    prof, young, flag = defaultdict(Counter), defaultdict(Counter), defaultdict(Counter)
    for q in cfg.codes:
        out = cfg.d("results", "combined") / f"{q}.tsv"
        ft = flag_track(cfg, q)
        with out.with_suffix(".tmp").open("w") as fh:
            fh.write("chrom\tstart\tend\tbranch\tsource\tyoungest\toldest\n")
            for c, x, y, b, src, lo, hi in combined_rows(cfg, q):
                fh.write(f"{c}\t{x}\t{y}\t{b}\t{src}\t{lo}\t{hi}\n")
                prof[q][(b, src)] += y - x
                if ft:
                    flag[q][b] += flagged_bp(ft, c, x, y)
                if b == "undated":
                    young[q][lo] += y - x
        out.with_suffix(".tmp").rename(out)
        print(f"[combine] {q}", flush=True)
    units = {q: unit_counts(cfg, q) for q in cfg.codes}
    nodes = sorted(tree.parent, key=lambda n: (-len(tree.leaves_of[n] & set(cfg.codes)), n))
    with (cfg.d("results") / "combined_branch.tsv").open("w") as fh:
        fh.write("branch\tn_genomes\tstepwise_mbp\ttile_fill_no_homolog_mbp\ttile_fill_in_range_mbp\tcopy_mbp\tgap_fill_mbp\t"
                 "combined_mbp\t"
                 "undated_youngest_here_mbp\tcombined_upper_mbp\t" +
                 "\t".join(f"{k}\t{k}_kbp_each" for k in UNIT_GAPS) + "\tflagged_mbp\tcombined_unflagged_mbp\n")
        for n in nodes:
            desc = sorted(tree.leaves_of[n] & set(cfg.codes))
            if not desc:
                continue
            k = len(desc) * 1e6
            a, b, c, cp, gf = (sum(prof[q][(n, s)] for q in desc) / k
                               for s in ("stepwise", "tile", "tile_in_range", "copy", "gap_fill"))
            u = sum(young[q][n] for q in desc) / k
            uc = []
            for kind in UNIT_GAPS:
                nu = sum(units[q][0][kind][n] for q in desc) / len(desc)
                ub = sum(units[q][1][kind][n] for q in desc) / len(desc)
                uc += [f"{nu:.1f}", f"{ub / nu / 1e3:.1f}" if nu else ""]
            fl = sum(flag[q][n] for q in desc) / k
            fh.write(f"{n}\t{len(desc)}\t{a:.2f}\t{b:.2f}\t{c:.2f}\t{cp:.2f}\t{gf:.2f}\t{a + b + c + cp + gf:.2f}\t{u:.2f}\t"
                     f"{a + b + c + cp + gf + u:.2f}\t"
                     + "\t".join(uc) + f"\t{fl:.2f}\t{a + b + c + cp + gf - fl:.2f}\n")
    with (cfg.d("results") / "combined_genome.tsv").open("w") as fh:
        fh.write("genome\tsd_mbp\tspecific_stepwise\tspecific_tile_fill\tspecific_combined\tshared_stepwise\t"
                 "shared_tile_fill\tshared_combined\tundated_range\tundated_no_homolog\tundated_tile_too_old\t"
                 "undated_pct\tspecific_flagged\tshared_flagged\tundated_flagged\n")
        for q in cfg.codes:
            p = prof[q]
            tot = sum(p.values())
            ss = p[(q, "stepwise")]
            st = p[(q, "tile")] + p[(q, "tile_in_range")] + p[(q, "copy")] + p[(q, "gap_fill")]
            hs = sum(v for (b, s), v in p.items() if s == "stepwise" and b != q)
            ht = sum(v for (b, s), v in p.items() if (s.startswith("tile") or s in ("copy", "gap_fill")) and b != q and b != "undated")
            ur, un, upar = p[("undated", "range")], p[("undated", "no_homolog")], p[("undated", "tile_too_old")]
            fs, fu = flag[q][q], flag[q]["undated"]
            fh.write("\t".join([q] + [C.mbp(v) for v in (tot, ss, st, ss + st, hs, ht, hs + ht, ur, un, upar)]
                               + [f"{100 * (ur + un + upar) / tot:.0f}"]
                               + [C.mbp(v) for v in (fs, sum(flag[q].values()) - fs - fu, fu)]) + "\n")
    print((cfg.d("results") / "combined_branch.tsv").read_text())


def per_sd(cfg):
    """One row per SD (merged SD interval, and SD unit), dated with the same rules as `annotate`."""
    out_dir = cfg.d("results", "per_sd")
    summary = []
    for q in cfg.codes:
        ann = Annotator(cfg, q)
        sd = C.sd_by_chrom(cfg, q)
        for kind, regions in (("sd_intervals", lambda ivs: [(a, b, b - a) for a, b in ivs]),
                              ("sd_units", lambda ivs: [(a, b, bp) for a, b, bp, _ in C.O.chain_iv(ivs, UNIT_GAPS["units_5kb"])])):
            n = Counter()
            with (out_dir / f"{q}.{kind}.tsv").open("w") as fh:
                fh.write("chrom\tstart\tend\tsd_bp\t" + ANN + "\n")
                for c in sorted(sd):
                    for a, b, sdbp in regions(sd[c]):
                        v = ann(c, a, b)
                        fh.write(f"{c}\t{a}\t{b}\t{sdbp}\t" + "\t".join(v) + "\n")
                        n[v[3]] += 1
                        if v[-1] and int(v[-1]) >= 50:
                            n[v[3] + "_flagged"] += 1
            summary.append((q, kind, n))
        print(f"[per-sd] {q}", flush=True)
    with (out_dir / "summary.tsv").open("w") as fh:
        fh.write("genome\tunit\tn\tspecies_specific\tshared\tundated\tspecies_specific_flagged\tshared_flagged\t"
                 "undated_flagged\n")
        for q, kind, n in summary:
            fh.write(f"{q}\t{kind}\t{n['species_specific'] + n['shared'] + n['undated'] + n['not_sd']}\t"
                     f"{n['species_specific']}\t{n['shared']}\t{n['undated']}\t{n['species_specific_flagged']}\t"
                     f"{n['shared_flagged']}\t{n['undated_flagged']}\n")
    print((out_dir / "summary.tsv").read_text())


# ------------------------------------------------------------------ annotate regions (e.g. SEDEF rows)
MIN_DATED = 0.5      # a region gets a main branch only when at least this share of its called bp is dated
ANN = ("bp\tcalled_bp\tdated_pct\tclass\tmain_branch\tmain_pct_of_dated\tbranch_genomes\tspecies_specific_pct\t"
       "shared_pct\tundated_pct\tstepwise_pct\ttile_pct\tcopy_pct\tgap_fill_pct\tsecond_branch\tbranches\tflagged_pct")
N_ANN = len(ANN.split("\t"))


class Annotator:
    """Dates any region of genome Q from the combined per-base calls (results/combined/<Q>.tsv).
    called_bp: bases in the pipeline's SD set (only those carry a call); % columns are of called_bp.
    class: 'undated' when < MIN_DATED of called_bp is dated (main_branch then blank; `branches` still lists the
    dated part), else species_specific or shared, whichever holds more bp; 'not_sd' without called bp.
    stepwise_pct / tile_pct / copy_pct: share of called bp dated by each method (tile = tile and tile_in_range;
    copy = per-copy gene-tree call, percopy.py; gap_fill = short undated run between identical calls).
    flagged_pct: share of called bp inside the genome's flag_regions (blank without them)."""

    def __init__(self, cfg, q):
        self.cfg, self.q = cfg, q
        self.tr = C.track([(r["chrom"], int(r["start"]), int(r["end"]), (r["branch"], r["source"])) for r in
                           csv.DictReader((cfg.d("results", "combined") / f"{q}.tsv").open(), delimiter="\t")])
        self.ft = flag_track(cfg, q)

    def __call__(self, c, a, b):
        per, src = Counter(), Counter()
        fl = 0
        for x, y, lab in C.walk(self.tr, c, a, b):
            if lab is not None:
                per[lab[0]] += y - x
                src[lab[1]] += y - x
                fl += flagged_bp(self.ft, c, x, y)
        called = sum(per.values())
        if not called:
            return [str(b - a), "0", "", "not_sd"] + [""] * (N_ANN - 4)
        und = per.pop("undated", 0)
        spec = per.get(self.q, 0)
        shared = called - und - spec
        dated = called - und
        top = per.most_common(2)
        if dated < MIN_DATED * called:
            cls, main = "undated", ""
        else:
            cls, main = ("species_specific" if spec >= shared else "shared"), top[0][0]
        pc = lambda v: f"{100 * v / called:.0f}"   # noqa: E731
        return [str(b - a), str(called), pc(dated), cls, main,
                f"{100 * top[0][1] / dated:.0f}" if main else "",
                str(len(self.cfg.tree.leaves_of[main] & set(self.cfg.codes))) if main else "",
                pc(spec), pc(shared), pc(und), pc(src["stepwise"]), pc(src["tile"] + src["tile_in_range"]), pc(src["copy"]),
                pc(src["gap_fill"]),
                (top[1][0] if len(top) > 1 else "") if main else "",
                ",".join(f"{k}:{v}" for k, v in per.most_common()), pc(fl) if self.ft else ""]


def annotate_file(cfg, q, path, out):
    """BED (chrom start end ...) or SEDEF/BEDPE (sides in columns 1-3 and 10-12): each region / side dated.
    Pairs also get same_main_branch (yes/no). Coordinates must be those of the HAL genome."""
    ann = Annotator(cfg, q)
    n = 0
    with open(path) as fi, open(out, "w") as fo:
        mode = None
        for line in fi:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if mode is None:
                mode = "pair" if len(f) >= 12 and f[10].isdigit() and f[11].isdigit() else "bed"
                if mode == "pair":
                    fo.write("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\t" +
                             "\t".join(f"side1_{x}" for x in ANN.split("\t")) + "\t" +
                             "\t".join(f"side2_{x}" for x in ANN.split("\t")) + "\tsame_main_branch\n")
                else:
                    fo.write("chrom\tstart\tend\t" + ANN + "\n")
            if mode == "pair":
                sides = []
                for c, a, b in ((f[0], f[1], f[2]), (f[9], f[10], f[11])):
                    sides.append(ann(c, int(a), int(b)) if c else [""] * N_ANN)
                same = "yes" if sides[0][4] and sides[0][4] == sides[1][4] else ("no" if sides[0][4] and sides[1][4] else "")
                fo.write("\t".join([f[0], f[1], f[2], f[9], f[10], f[11]] + sides[0] + sides[1] + [same]) + "\n")
            else:
                fo.write("\t".join([f[0], f[1], f[2]] + ann(f[0], int(f[1]), int(f[2]))) + "\n")
            n += 1
    return n


def annotate(cfg, args):
    if args[0] == "all":
        for q in cfg.codes:
            out = cfg.d("results", "annotated_sedef") / f"{q}.tsv"
            n = annotate_file(cfg, q, cfg.species[q]["sedef"], out)
            print(f"[annotate] {q}: {n} SEDEF rows -> {out}", flush=True)
        return
    q, target = args[0], args[1]
    if q not in cfg.codes:
        sys.exit(f"unknown genome {q}; genomes: {' '.join(cfg.codes)}")
    if ":" in target and not Path(target).exists():
        c, r = target.rsplit(":", 1)
        a, b = (int(x.replace(",", "")) for x in r.split("-"))
        vals = Annotator(cfg, q)(c, a, b)
        for k, v in zip(["region"] + ANN.split("\t"), [target] + vals):
            print(f"{k:20} {v}")
        return
    out = args[2] if len(args) > 2 else str(Path(target).with_suffix("")) + f".{q}.sdplace.tsv"
    n = annotate_file(cfg, q, target, out)
    print(f"[annotate] {q}: {n} rows -> {out}")


def validate(cfg, regions):
    """Regions of known duplication age: % of dated SD bp on the expected branch, older, younger; % undated."""
    by_q = defaultdict(list)
    with open(regions) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            by_q[r["genome"]].append(r)
    with (cfg.d("results") / "validation.tsv").open("w") as fh:
        fh.write("name\tgenome\texpected\tsd_bp\tmethod\tundated_pct\texpected_pct\tolder_pct\tyounger_pct\ttop_call\n")
        for q, regs in by_q.items():
            idx = {n: i for i, n in enumerate(cfg.path_to_root(q))}
            sd = C.sd_by_chrom(cfg, q)
            comb = C.track([(c, x, y, b) for c, x, y, b, *_ in combined_rows(cfg, q)])
            step = stepwise_track(cfg, q)
            tracks = {"combined": comb, "stepwise": step, "tiles": tile_track(cfg, q)}
            for r in regs:
                exp = r["expected_branch"]
                if exp not in idx:
                    sys.exit(f"{r['name']}: '{exp}' is not on the path from {q} to the root")
                pcs = C.sd_pieces(sd, r["chrom"], int(r["start"]), int(r["end"]))
                sdbp = sum(b - a for a, b in pcs)
                for nm, tr in tracks.items():
                    cl, calls = Counter(), Counter()
                    for a, b in pcs:
                        for x, y, l in C.walk(tr, r["chrom"], a, b):
                            if l not in idx:
                                cl["undated"] += y - x
                            else:
                                cl["expected" if l == exp else "older" if idx[l] > idx[exp] else "younger"] += y - x
                                calls[l] += y - x
                    d = sdbp - cl["undated"]
                    pc = lambda v, t: f"{100 * v / t:.0f}" if t else ""   # noqa: E731
                    fh.write(f"{r['name']}\t{q}\t{exp}\t{sdbp}\t{nm}\t{pc(cl['undated'], sdbp)}\t{pc(cl['expected'], d)}\t"
                             f"{pc(cl['older'], d)}\t{pc(cl['younger'], d)}\t{calls.most_common(1)[0][0] if calls else ''}\n")
    print((cfg.d("results") / "validation.tsv").read_text())


def main():
    a = sys.argv[1:]
    cfgp = None
    if a and a[0] == "-c":
        cfgp, a = a[1], a[2:]
    cfg = C.load(cfgp)
    if a and a[0] == "validate":
        validate(cfg, a[1])
    elif a and a[0] == "per-sd":
        per_sd(cfg)
    elif a and a[0] == "annotate":
        annotate(cfg, a[1:])
    else:
        combine(cfg)


if __name__ == "__main__":
    main()
