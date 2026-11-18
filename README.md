# sdplace

Dates segmental duplications (SDs) on a species tree: for every SD base of every genome, the branch on
which that sequence became SD.

Inputs
- A Cactus HAL containing every genome.
- A soft-masked genome FASTA of each genome
- One SEDEF SD file per genome, filtered (identity > 90%, length > 1 kbp, satellite < 70%).
- A rooted species tree in Newick with every internal node named

Software on the PATH: `python3` (>= 3.9) with `pysam` and `mappy`; `samtools`; `minimap2`; `blastn` (BLAST+);
`halLiftover` and `halStats`.

## Prepare your config file

Copy `config.example.json` and edit it. Relative paths are relative to the config file.

| Key | Meaning |
|---|---|
| `hal`, `outdir`, `tree` | The alignment, output directory, the species tree |
| `hal_command_prefix` | Command prefix for the HAL tools |
| `species.<CODE>` | `genome` (name in the HAL), `fasta`, `sedef`, optional `chrom_map`, `exclude_chroms` (e.g. chrY) |
| `minimap2_anchor_clades` | Tree nodes inside which unlifted rescue anchors may be placed with minimap2. Use nodes whose genomes are within about 10 My of each other (sister species, a genus). Default: none (Cactus only) |
| `asm20_clades` | Tree nodes inside which SD is aligned with minimap2 asm20 (up to roughly 20% divergence); other pairs use masked blastn. Default: minimap2 for every pair |
| `blast_priority` | Optional order in which outgroup genomes are tried in the tile blastn stage (default: config order) |
| `tile_age_cap` | Optional: `node_ages_my` (age of every internal node), `substitution_rate_per_site_per_my` (default 0.001), `margin` (default 2). A tile call is set aside (undated, `tile_too_old`) when its branch is entirely older than margin x (1 - identity) / (2 x rate) My, using the base's most diverged SEDEF paralog. Without it, tile calls are always used |
| `sedef_identity_column` | Column of the SEDEF files holding pair identity (default 24, SEDEF `fracMatch`) |
| `threads` | CPUs for the multi-threaded tasks |
| `slurm` | `partition`, `max_parallel`, `extra_args` (list of extra `sbatch` arguments) |

## Running the pipeline

```bash
python3 sdplace.py -c config.json status          # get status
python3 sdplace.py -c config.json run --slurm     # everything on slurm arrays
python3 sdplace.py -c config.json run --jobs 8    # everything on local
python3 sdplace.py -c config.json run align place combine   # specific stages, see below
python3 sdplace.py -c config.json run --dry-run
```

Each stage skips outputs that exist, so rerunning resumes. To redo a stage, delete its outputs.
With `--slurm` the driver waits for each array, so you might want to submit it as a job.

| Stage | What it does | Output (under `outdir`) |
|---|---|---|
| `prepare` | Merges SD; units (SD within 5 kbp), loci (SD within 100 kbp) with 100 kbp flanks; 1 kbp tiles | `prep/`, `tiles/bed/` |
| `lifts` | Cactus lifts (`halLiftover --noDupes`) of flanks, unit anchors and tiles, one task per genome pair and kind | `lifts/`, `tiles/lifts/` |
| `tile_cactus` | Per tile and other genome: SD there, sequence there but not SD, locus there without it, or cannot tell | `tiles/cactus/` |
| `tile_queries`, `tile_map`, `tile_mm2calls` | minimap2 asm20 for tiles Cactus cannot call | `tiles/queries/`, `tiles/paf/`, `tiles/mm2/` |
| `tile_blast` | Masked blastn inside the flank-anchored window for tiles still unresolved; final tile calls | `tiles/classified/` |
| `homologs` | Homolog of each locus = interval between its lifted flanks | `stepwise/homologs_locus.tsv` |
| `rescue_anchors`, `rescue_lift`, `rescue_mm2`, `rescue_homologs` | Loci whose flanks do not lift (satellite): outer windows walking out to 1 Mbp and interior non-SD gaps, lifted with Cactus, minimap2 only inside `minimap2_anchor_clades` | `stepwise/rescue/`, `stepwise/homologs_rescue.tsv` |
| `align` | Each genome's SD aligned to the SD inside each homolog (minimap2 asm20 or masked blastn) | `stepwise/aln/` |
| `place` | Stepwise call per base, walking from the genome to the root | `stepwise/placed/`, `stepwise/stepwise_genome.tsv` |
| `combine` | Stepwise call where resolved, tile call otherwise | `results/` |
| `per_sd` | One row per SD (merged interval and unit): species-specific, shared or undated | `results/per_sd/` |

## Results

- `results/combined_branch.tsv`: one row per branch with the SD Mbp per genome, SD units, averaged over descendant genomes. This is the file you need if you want to just look at how much Mbp of SD (place-able) sequence originated on a branch.
- `results/combined_genome.tsv`: per genome, species-specific, shared and undated Mbp. 
- `results/per_sd/<CODE>.sd_intervals.tsv` and `<CODE>.sd_units.tsv` (stage `per_sd`): one row per SD, either
  each merged SD interval or each SD unit (SD within 5 kbp). `class` is `species_specific`, `shared` or `undated`,
  whichever holds most of the SD's bp; `main_branch` is the branch holding most of its dated bp and
  `branch_genomes` the number of genomes under it; the same columns and rules as `annotate`. This is probably the file you need if you want to look at which branches different SDs originated from.
- `results/per_sd/summary.tsv` counts SDs per class and genome.


After `run`, any region of a genome can be dated from the per-base calls (coordinates of the HAL genome):
```bash
python3 sdplace.py -c config.json annotate HUMAN chr1:405368-408254  
python3 sdplace.py -c config.json annotate HUMAN regions.bed [out.tsv]                 
```

## Limitations

- Results depend on the Cactus alignment and on which genomes are present: a branch can only be told from
  its parent when an outgroup genome has the locus.
- SD flanked by satellite or other unalignable repeats often has no matched locus outside its sister
  species; it is then dated by the tile call or left undated (13-39% of a genome's SD in a set of primates).
- Tile-filled calls lean old where Cactus aligns to a paralog.
- Sequence comparison across very deep splits is unreliable (minimap2 asm20 is made for up to ~20%
  divergence); the deepest branches are floors.
- These scripts were put into a pipeline using AI agents. If you face any issues, please let me know!