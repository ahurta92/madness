# Excited-state benchmark suite (legacy vs v1)

Times each codebase on a ladder of systems, all running the same 4-root
TDA + RPA at protocol [1e-6], starting from a freshly-converged
moldft archive. The point is to measure the **rotate-vs-recompute gap**:
v1 carries the gamma matrix through subspace rotations, legacy recomputes
it. The gap should grow with `n_occ × n_states`, so we ramp `n_occ`.

## Layout

```
benchmarks/
├── README.md              (this file)
├── h2o/                   (n_occ=5, ~2 hr / node)
│   ├── benchmark.slurm
│   ├── response_tda.in
│   └── response_rpa.in
├── c2h4/                  (n_occ=8, ~6 hr / 2 nodes)
│   ├── benchmark.slurm
│   ├── response_tda.in
│   └── response_rpa.in
└── c6h6/                  (n_occ=21, ~24 hr / 4 nodes @ 4 tasks/node)
    ├── benchmark.slurm
    ├── response_tda.in
    └── response_rpa.in
```

The ground-state geometries and `moldft.in` files come from the existing
fixtures in `../systems/{h2o_hf, c2h4_hf, c6h6_hf}/`.

## What each `benchmark.slurm` does

1. Stages `moldft.in` into a per-job work dir under `$WORK_BASE`.
2. Runs **moldft** once to converge the ground state at protocol [1e-4, 1e-6].
3. For each `(code, mode) ∈ {legacy, v1} × {tda, rpa}`:
   - Symlinks the moldft archive into `run_<code>_<mode>/`.
   - Copies `response_<mode>.in` (named `input` for legacy, `response.in`
     for v1).
   - Runs the binary, captures wall-clock + iteration count + final omegas
     into `timing.csv`.
4. Prints a side-by-side summary at the end of the SLURM stdout.

## Before submitting

The SLURM scripts have **two edits** you almost certainly need:

```bash
# top of each benchmark.slurm:
#SBATCH --partition=long-40core    # <-- your seawulf partition
#SBATCH --account=stark            # <-- your account
```

You may also need to uncomment the `module load` lines if your build
needs them. The defaults assume:
- `$REPO=$HOME/Projects/madness`  (overridable via env)
- `$BUILD=$REPO/build`             (overridable via env)
- `$WORK_BASE=/gpfs/scratch/$USER/madness_es_bench` (overridable)

## Build prerequisites

On Seawulf, before submitting:

```bash
cd $REPO/build
make moldft -j8
make molresponse_legacy -j8
make molresponse        -j8           # this is v1
```

## Submit

```bash
cd src/apps/molresponse_v3/tests/fixtures/benchmarks/h2o   && sbatch benchmark.slurm
cd ../c2h4                                                 && sbatch benchmark.slurm
cd ../c6h6                                                 && sbatch benchmark.slurm
```

Each job is independent. C6H6 will take the longest; consider running
H2O first as a smoke test before committing C6H6's 24-hour reservation.

## After they finish

Each job leaves a `timing.csv` in its work dir with this schema:

```
system,code,mode,wall_seconds,iters,omega1,omega2,omega3,omega4,exit
```

A simple comparison across systems:

```bash
WORK_BASE=/gpfs/scratch/$USER/madness_es_bench
( head -1 $WORK_BASE/h2o_*/timing.csv | head -1
  for sys in h2o c2h4 c6h6; do
    tail -n+2 $WORK_BASE/${sys}_*/timing.csv 2>/dev/null
  done ) | column -t -s,
```

You should see `wall_seconds` for `code=v1` consistently below `code=legacy`
within the same `(system, mode)` pair. The ratio is the rotate-vs-recompute
speedup — the headline number we want for the design discussion.

## What we're looking for

The rotation savings live inside `compute_gamma_full` (the J/K/XC build).
With protocol [1e-6] and 4 roots, expected costs scale roughly as:

| n_occ | gamma cost  | rotate-savings expectation |
|-------|-------------|-----------------------------|
| 5  (H2O)  | small        | ~5-15% (small but measurable) |
| 8  (C2H4) | medium       | ~15-30%                       |
| 21 (C6H6) | dominates    | ~40-60% (this is the target)  |

If the gap on C6H6 is < 30% something is off and we re-think before
porting the rotate-everything pattern to v3. If it's > 50%, port it.
