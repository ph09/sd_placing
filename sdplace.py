#!/usr/bin/env python3
"""Driver
"""
import argparse
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402


def py(script, cfg, *args):
    return " ".join(["python3", shlex.quote(str(HERE / script)), "-c", shlex.quote(str(cfg.path))] + [str(a) for a in args])


def stages(cfg):
    """[(name, kind, cpus, mem_gb, hours, todo())]: todo -> [command] still to run. kind: serial | array"""
    cs = cfg.codes
    pairs = [(q, t) for q in cs for t in cs if q != t]
    out = cfg.out
    th = cfg.threads

    def missing(paths):
        return [p for p in paths if not Path(p).is_file()]

    def partners(t):
        return [q for q in cs if q != t and cfg.share_clade(q, t, cfg.mm2_anchor_clades)]

    return [
        ("prepare", "serial", 2, 16, 2, lambda: (
            [py("prepare.py", cfg)] if missing(cfg.prep(q, "loci.tsv") for q in cs) else []) + (
            [py("tiles.py", cfg, "prepare")] if missing(out / "tiles" / "bed" / f"{q}.bed" for q in cs) else [])),
        ("lifts", "array", 2, 60, 6, lambda: (
            [py("lifts.py", cfg, k, q, t) for k in ("flanks", "anchors") for q, t in pairs
             if not cfg.lift(k, q, t).is_file()] +
            [py("tiles.py", cfg, "lift", q, t) for q, t in pairs
             if not (out / "tiles" / "lifts" / f"{q}__TO__{t}.tsv").is_file()])),
        ("tile_cactus", "array", 2, 48, 4, lambda: [
            py("tiles.py", cfg, "cactus", q) for q in cs if not (out / "tiles" / "cactus" / f"{q}.tsv").is_file()]),
        ("tile_queries", "serial", 2, 16, 4, lambda: (
            [py("tiles.py", cfg, "queries")] if missing(out / "tiles" / "queries" / f"{q}.fa" for q in cs) else [])),
        ("tile_map", "array", th, 96, 12, lambda: [
            py("tiles.py", cfg, "map", t, th) for t in cs if not (out / "tiles" / "paf" / f"to_{t}.paf").is_file()]),
        ("tile_mm2calls", "serial", 2, 96, 4, lambda: (
            [py("tiles.py", cfg, "mm2calls")] if missing(out / "tiles" / "mm2" / f"{q}.tsv" for q in cs) else [])),
        ("tile_blast", "array", th, 96, 12, lambda: [
            py("tiles.py", cfg, "blast", q, th) for q in cs
            if not (out / "tiles" / "classified" / f"{q}.tsv").is_file()]),
        ("homologs", "serial", 2, 32, 2, lambda: (
            [py("stepwise.py", cfg, "homologs")] if missing([out / "stepwise" / "homologs_locus.tsv"]) else [])),
        ("rescue_anchors", "serial", 2, 32, 2, lambda: (
            [py("stepwise.py", cfg, "anchors")]
            if missing(out / "stepwise" / "rescue" / "anchors" / f"{q}__TO__{t}.bed" for q, t in pairs) else [])),
        ("rescue_lift", "array", th, 100, 12, lambda: [
            py("stepwise.py", cfg, "rescue-lift", q, th) for q in cs
            if missing(out / "stepwise" / "rescue" / "cactus" / f"{q}__TO__{t}.tsv" for t in cs if t != q)]),
        ("rescue_mm2", "array", th, 100, 12, lambda: [
            py("stepwise.py", cfg, "rescue-mm2", t, th) for t in cs
            if missing(out / "stepwise" / "rescue" / "mm2" / f"{q}__TO__{t}.tsv" for q in partners(t))]),
        ("rescue_homologs", "serial", 2, 80, 4, lambda: (
            [py("stepwise.py", cfg, "rescue-homologs")] if missing([out / "stepwise" / "homologs_rescue.tsv"]) else [])),
        ("align", "array", th, 100, 12, lambda: [
            py("stepwise.py", cfg, "align", q, th) for q in cs if not (out / "stepwise" / "aln" / f"{q}.tsv").is_file()]),
        ("place", "serial", 2, 80, 2, lambda: (
            [py("stepwise.py", cfg, "place")] if missing(out / "stepwise" / "placed" / f"{q}.tsv" for q in cs) else [])),
        ("percopy_select", "serial", 2, 80, 2, lambda: (
            [py("percopy.py", cfg, "select")] if missing(out / "percopy" / "queries" / f"{q}.fa" for q in cs) else [])),
        ("percopy_map", "array", th, 120, 12, lambda: [
            py("percopy.py", cfg, "map", t, th) for t in cs if not (out / "percopy" / "paf" / f"to_{t}.paf").is_file()]),
        ("percopy_trees", "array", 32, 96, 12, lambda: [
            py("percopy.py", cfg, "trees", q, 32) for q in cs
            if not (out / "percopy" / "calls" / f"{q}.tsv").is_file()]),
        ("combine", "serial", 2, 80, 2, lambda: (
            [py("combine.py", cfg)] if missing([out / "results" / "combined_branch.tsv"]) else [])),
        ("per_sd", "serial", 2, 32, 2, lambda: (
            [py("combine.py", cfg, "per-sd")] if missing([out / "results" / "per_sd" / "summary.tsv"]) else [])),
    ]


def run_local(cmds, jobs, log_dir, name):
    def one(ic):
        i, c = ic
        with (log_dir / f"{name}_{i}.log").open("w") as log:
            return subprocess.call(c, shell=True, stdout=log, stderr=subprocess.STDOUT)
    with ThreadPoolExecutor(max(1, jobs)) as ex:
        return sum(1 for rc in ex.map(one, enumerate(cmds)) if rc)


def run_slurm(cfg, cmds, name, cpus, mem, hours, log_dir):
    f = log_dir / f"{name}.cmds"
    f.write_text("\n".join(cmds) + "\n")
    s = cfg.slurm
    cmd = ["sbatch", "--wait", f"--job-name=sdplace_{name}", f"--cpus-per-task={cpus}", f"--mem={mem}G",
           f"--time={int(hours)}:00:00", f"--array=0-{len(cmds) - 1}%{int(s.get('max_parallel', 60))}",
           f"--output={log_dir}/{name}_%A_%a.out"]
    if s.get("partition"):
        cmd.append(f"--partition={s['partition']}")
    cmd += list(s.get("extra_args", []))
    cmd.append(f"--wrap=eval \"$(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" {shlex.quote(str(f))})\"")
    return subprocess.call(cmd)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", required=True)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("status")
    p = sub.add_parser("run")
    p.add_argument("stages", nargs="*", help="default: all, in order")
    p.add_argument("--slurm", action="store_true")
    p.add_argument("--jobs", type=int, default=1, help="local: tasks in parallel")
    p.add_argument("--dry-run", action="store_true", help="print the commands that would run")
    p = sub.add_parser("validate")
    p.add_argument("regions")
    p = sub.add_parser("annotate", help="date regions from the combined calls (after `run`)")
    p.add_argument("what", nargs="+", help="all | GENOME chrom:start-end | GENOME FILE [OUT]")
    args = ap.parse_args()
    cfg = C.load(args.config)
    st = stages(cfg)
    if args.cmd == "status":
        for name, kind, *_, todo in st:
            n = len(todo())
            print(f"{name:16} {'done' if not n else f'{n} task(s) to run'}")
        return
    if args.cmd == "validate":
        sys.exit(subprocess.call(py("combine.py", cfg, "validate", args.regions), shell=True))
    if args.cmd == "annotate":
        sys.exit(subprocess.call(py("combine.py", cfg, "annotate", *[shlex.quote(x) for x in args.what]), shell=True))
    names = [s[0] for s in st]
    for x in args.stages:
        if x not in names:
            sys.exit(f"unknown stage '{x}'; stages: {' '.join(names)}")
    log_dir = cfg.d("logs", "sdplace")
    for name, kind, cpus, mem, hours, todo in st:
        if args.stages and name not in args.stages:
            continue
        for attempt in (1, 2):
            cmds = todo()
            if not cmds:
                break
            print(f"[{name}] {len(cmds)} task(s)" + (" (retry)" if attempt == 2 else ""), flush=True)
            if args.dry_run:
                print("\n".join(cmds))
                break
            if args.slurm and kind == "array":
                run_slurm(cfg, cmds, name, cpus, mem, hours, log_dir)
            else:
                run_local(cmds, args.jobs if kind == "array" else 1, log_dir, name)
        if not args.dry_run and todo():
            sys.exit(f"[{name}] did not finish: see {log_dir}/{name}_*")
        if not args.dry_run:
            print(f"[{name}] done", flush=True)


if __name__ == "__main__":
    main()