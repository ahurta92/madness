# MRA operation benchmarks

A microbenchmark harness for measuring how each core MADNESS MRA operation scales
with the **number of functions `n`** and the **protocol `(thresh, k)`**, plus MPI
rank count `P` and threads. It is the measurement foundation for the performance
models in [`docs/parallel_runtime_guide/`](../../../docs/parallel_runtime_guide/)
and the I/O / exchange work that builds on them.

Files: `bench_mra_ops.cc` (driver), `analyze.py` (summarizer), this README.

---

## The methodology (what we measure and why)

We treat each operation as an isolated experiment and sweep a grid:

```
axes:   operation  x  n (#functions)  x  k  x  thresh   [ x P x threads ]
hold fixed within a point: cell L, function content (seeded generator),
                           operator (built once per protocol)
```

**Controlled inputs.** Functions are sums of a few Gaussians at *deterministic,
seeded* centers (`GaussianBlob`), with `special_points()` forcing refinement at
the centers (orbital-cusp-like). Because the seed depends only on the function
index, the tree for a given `n` is reproducible across protocols — so the only
things changing `N_leaf` are `k` and `thresh`, never run-to-run noise.

**What each record captures** (one JSONL line per `(op, n, k, thresh, P, threads)`):

| Field | Meaning | Model term it feeds |
|-------|---------|---------------------|
| `time_s.median/min/max` | wall time over `reps` (warmup discarded) | the thing we model |
| `tree.nodes` / `tree.leaves` | global node / leaf count of the operand set | problem size `N`, `N_leaf` |
| `tree.nodes_rank_min/max`, `imbalance` | per-rank local node spread, `max·P/total` | load-balance factor `φ` |
| `tree.max_depth` | deepest level | `L` |
| `coeff_bytes_est` | `leaves · k³ · 8` | memory model (companion §11) |
| `rmi.nmsg_sent/nbyte_sent/...` | global RMI traffic during the op | `T_comm = n_msg·α + bytes·β` |

These map one-to-one onto `T_op = T_compute + T_comm + T_sync` and the memory
model. The point is not just timing — it is timing **alongside** the counts that
explain the timing, so you can fit and *predict* rather than re-measure.

**Timing discipline** (in `time_op`):
- One untimed **warmup** run (first-touch allocation, operator caches).
- Inputs are (re)built by `prepare()` *outside* the timed region; only the op
  call is between the fences.
- Each timed run is `fence(); t0; op; fence(); t1` so the time includes the
  operation's own quiescence (its remote tasks completing) — the realistic cost.
- RMI counters are sampled inside the fences and **summed across ranks**, so the
  reported traffic is the global total attributable to the op. (The closing
  fence's own `O(log P)` control messages are included; they are negligible next
  to operation traffic and constant across ops.)
- We report the **median** of `reps` runs.

**Operations covered** (vector form — the realistic molresponse path):
`project` (construct from functor), `compress`, `reconstruct`, `truncate`,
`gaxpy` (`A += B`), `multiply` (`A[i]·B[i]`), `matrix_inner` (`n×n` Gram), and
`apply` (BSH or Coulomb convolution — the heavy op). Each is fed inputs in the
tree state it expects (e.g. `apply`/`multiply` get reconstructed inputs;
`gaxpy`/`matrix_inner` get compressed) so the timed region is the op's core work,
not a hidden state conversion.

---

## Build

```bash
# from a configured build dir (see CLAUDE.md for the cmake line)
ninja bench_mra_ops
```

`src/apps/bench` is wired into `src/apps/CMakeLists.txt`. The target links
`MADchem` (for the operator constructors) and installs to the bin dir.

## Run

```bash
# single rank, sweep protocol and n
./bench_mra_ops --k=6,8,10 --thresh=1e-4,1e-6 --n=1,4,16 --reps=3 --out=bench.jsonl

# distributed: 4 ranks x 8 threads each
MAD_NUM_THREADS=8 mpiexec -np 4 ./bench_mra_ops \
    --k=6,8,10 --thresh=1e-4,1e-6 --n=1,4,16,32 --reps=3 \
    --operator=bsh --L=20 --out=bench_p4.jsonl
```

Options (all have defaults): `--k`, `--thresh`, `--n` (comma lists); `--L`,
`--reps`, `--gaussians`, `--width`; `--operator=bsh|coulomb`; `--out`.

## Analyze

```bash
python3 analyze.py bench.jsonl            # all ops
python3 analyze.py bench.jsonl --op=apply # one op
```

`analyze.py` reports, per operation, the fitted exponents `time ~ n^p` (at fixed
`k,thresh`) and `time ~ k^q` (at fixed `n,thresh`), plus the communication and
imbalance per point. Expected sanity checks against the models in the guide:

- `compress`/`reconstruct`/`multiply`/`truncate`: `p ≈ 1` in `n`, and `q`
  approaching the transform exponent (toward `~k^{d+1}` per leaf, modulated by how
  `N_leaf` grows with `k`).
- `matrix_inner`: `p ≈ 2` in `n` (it is the `n×n` Gram).
- `apply`: the steepest `q`; watch `rmi.nbyte_sent` and `imbalance` climb with `P`.

## How to use it for the three workstreams

1. **Calibrate the performance model.** Run on H2O/CH3OH-scale inputs, fit `R, α,
   β, φ` from `time` vs the measured counts, then predict C6H6/naphthalene before
   launching (compare to the `mul_sparse` table in `CLAUDE.md`).
2. **HDF5 I/O.** Add `save`/`load` timing points (same harness shape) to measure
   the single-writer ceiling and, later, the parallel-HDF5 speedup.
3. **Exchange (K).** Add a `K`-specific record (pairwise `multiply` + `apply`)
   once the refactor lands, reusing the same metrics so it is directly comparable.

> Note: this harness has not yet been compiled in this workspace (no build dir was
> present). Build it in a configured tree and report any API drift.
