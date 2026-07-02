# 29 — Performance model: instrumentation + cost model (perf-model thread)

> **Thread:** `perf-model` (branch `perf-model`, off `molresponse-feature-next`).
> The *measurement arm* of the release effort. Mandate, test recipe, contracts:
> the **perf-model** brief in `madness_studies/RELEASE_STATUS.md`. This doc is the
> thread's design anchor; the running log is the `## perf-model log` section of
> `docs/00_status.md`.

## Mandate (what this thread delivers)

1. **Per-phase timers/counters in the core** — `apply` / `compress` /
   `reconstruct` / `multiply` / `inner` / `gaxpy` / `truncate`, plus the
   response-level phases **exchange** (γ build) and **projection** (`Q·v`).
2. **A machine-readable profile** emitted per run (JSON), with a pinned schema.
3. **A cost model** that predicts wall-time from problem shape
   (n_occupied, k, box, ranks), fit on a small sweep and checked
   predicted-vs-measured.

**Exposes (inter-thread contract):** the meter API + the profile **schema** that
`exchange` reports its Tx/tile counts + phase timings into, and that
`parallel-runtime` uses as the quantitative tiebreaker for the doc-24 (persistent
subworlds) vs doc-25 (pmaps/`Group`/`LoadBalanceDeux`) fork. The schema is pinned
in `docs/operator_contracts.md` (§ "Performance profile schema").

**Do-not-touch / hard contracts:**
- ❌ Numerics / solver convergence. Instrumentation must be **zero-effect when
  off** (compile flag *and* runtime env gate; see below).
- ❌ Exchange's reference kernels `compute_V0x` / `compute_gamma` / `compute_E0x`
  — they are exchange-thread's gate-0 oracle. We may wrap them in named meter
  blocks but must not change their numerics or control flow.
- Propose a concrete diff + get approval before non-trivial core/solver edits;
  the USER runs builds/solves on the allocation.

## Design decision: reuse core `WorldProfile` (not a new meter)

Core MADNESS already ships the meter we need: `WorldProfile`
(`src/madness/world/worldprofile.{h,cc}`).

- **Already gated, already zero-cost when off.** Compile flag
  `ENABLE_WORLD_PROFILE` → `WORLD_PROFILE_ENABLE` (CMake, default **OFF**). When
  off, the `PROFILE_FUNC` / `PROFILE_BLOCK` / `PROFILE_MEMBER_FUNC` macros expand
  to nothing — literally no code. This *is* the "zero-effect when off" contract.
- **Already instrumented.** The core ops carry `PROFILE_*` macros today:
  `mraimpl.h` (38), `funcimpl.h` (25), `operator.h` (13), `vmra.h` (39) — covering
  apply / compress / reconstruct / project / truncate.
- **Already captures §12's quantities.** Per call-site, parallel-reduced
  (sum / min / max + which rank): exclusive & inclusive **CPU time**, **call
  count**, and **messages + bytes sent/received** (`WorldProfileEntry`,
  `worldprofile.h:84-125`; per-thread accumulation in `WorldProfileObj`,
  `worldprofile.cc:407-461`). These are exactly the inputs the §9 cost model
  needs for `T_compute` (cpu + count), `T_comm` (msgs + bytes), and the imbalance
  factor φ (max vs sum across ranks).
- `WorldProfile::print(world)` (`worldprofile.cc:295`) already does the binary-tree
  reduction to rank 0; it is invoked from `World::print_stats` →
  `madqc.cpp:227`. We mirror that reduction for a JSON dump.

A v3-local wall-timer was rejected: it would duplicate this facility and could not
recover the intra-core msg/byte attribution (the `T_comm` inputs) without
re-instrumenting the core anyway. **Wall** time is recovered by *joining* the fine
per-phase CPU/comm profile with the existing coarse wall layer
(`StateMetrics.wall_s` per state/protocol; `PROTOCOL_START/DONE` lines) — see PM-3.

## Increment plan

### PM-1 — JSON emitter + pinned schema  *(first; smallest, additive, no numerics)*
- Add `WorldProfile::dump_json(World&, const std::string& path)` to
  `worldprofile.{h,cc}`. It performs the **same** tree reduction as `print`, but
  on rank 0 serializes the reduced `std::vector<WorldProfileEntry>` to JSON
  instead of formatting a human table. Refactor the reduction out of `print` into
  a private `reduce_to_root(World&) -> std::vector<WorldProfileEntry>` so `print`
  and `dump_json` share one gather and neither double-reduces the live `items`.
- **Double gate:** body compiled only under `WORLD_PROFILE_ENABLE`; at the call
  site, emit only if env `MADQC_PROFILE_JSON` is set (→ its value is the output
  path, else default `<cwd>/perf_profile.json`). So a profile-enabled build is
  still silent unless explicitly asked.
- **Call site:** alongside the existing `print_stats(world)` at madqc shutdown
  (`madqc.cpp:227`) and at the `run_response` seam end, so both the madqc and the
  standalone `run_response` app paths emit.
- **Schema:** pinned in `operator_contracts.md`. Shape-agnostic in the core
  emitter (WorldProfile knows nothing about orbitals); an optional `context`
  object is filled by the v3 caller with problem shape (molecule, n_occ, k,
  thresh, protocol, P) so PM-3 can join without a separate file.

### PM-2 — canonical phase taxonomy + response-level meters
- Map the raw `__FUNCTION__`-keyed entries onto canonical phases
  `{apply, compress, reconstruct, multiply, inner, gaxpy, truncate, exchange,
  projection}` via a small name→phase table in the fit script (keeps the core
  edit zero). Unmapped entries pass through as `phase=other`.
- Add two **named** `PROFILE_BLOCK`s at the v3 response phases that are not core
  ops: `rs_exchange_gamma` around `compute_gamma` (`kernels/full.hpp`) and
  `rs_projection` around `rs::project` (`response_space_ops.hpp`). **This is the
  meter API `exchange` reports into** — the exchange thread adds its Tx/tile
  counters inside the same blocks so they aggregate under `exchange`.

### PM-3 — cost-model fit
- Offline script `madness_studies/refs/perf_model_fit.py` (sibling of
  `study_analyze.py`). Ingests `perf_profile.json` across a small (k, n_occ, P)
  sweep, joins per-run problem shape (from the `context` block or
  `response_metadata.json` in the same dir), and fits §9's
  `T = T_compute + T_comm + T_sync`:
  - `T_compute ≈ (N_busiest · c_node) / (R · n_threads_eff)`, `c_node ∝ k^{d+1}`
  - `T_comm   = n_msg · α + V_bytes · β`
  - imbalance `φ = N_max_rank · P / N`
  Solve for machine constants `R, α, β` and report **predicted vs measured wall**
  for held-out shapes (target: predict C6H6 / naphthalene before launch — the
  `mul_sparse` table in the top-level `CLAUDE.md`).

## Profile schema (v1) — see `operator_contracts.md` for the pinned copy
Per run, rank-0 JSON. Each per-phase stat is `{sum, min, max, pmin, pmax}` faithful
to `WorldProfileEntry`. Run-level: `world_size`, `total_cpu_s`, `total_wall_s`,
`overhead_s_per_call`. Optional `context` (problem shape) filled by the caller.

## Test recipe (the brief's gate)
```
cm_use perf-model
cm_build           # MUST add -DENABLE_WORLD_PROFILE=ON for this thread
MADQC_PROFILE_JSON=$PWD/perf_profile.json cm_run h2o
# → confirm perf_profile.json emits; numbers stable across 2 runs
# → fit on a small (k, n_occ) sweep; check predicted vs measured wall
```
Core-lib edits (worldprofile.cc/.h) need a **full `ninja`**, not just cm targets
(they affect all of MADNESS).

## Status
PM-1 design complete; concrete diff proposed, awaiting approval + cm_build run.
PM-2 / PM-3 follow once PM-1's profile emits and is stable.
