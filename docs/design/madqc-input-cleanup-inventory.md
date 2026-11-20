# Remaining response inputs after this cleanup

Retired public inputs (rejected and omitted from listings): accuracy.*, convergence.use_default_criteria, progress.*, diagnostic.iterations_per_rung, stall.*, seed.start_rung, seed.es_warmup, excited.maxiter, excited.guess_max_iter, dconv. Inactive internal fields remain for private legacy kernel harnesses. Ordinary madqc response uses the shared controller.

| Remaining group | Why still present / next decision |
| --- | --- |
| convergence.* | Shared criteria/tolerances for all response states; provisional absolute defaults are documented in the response README. |
| excited.guess, excited.guess_basis, excited.maxsub | Retained guess selection and KAIN controls. Cold guesses select from the existing larger Ritz trial space; no BSH/KAIN warmup updates. The old guess iteration budget is retired. |
| kain.min_residual | Existing residual-dependent KAIN gate used by older solver paths; shared FD advancement currently uses its separate KAIN adapter. Decide its shared applicability explicitly. |
| engine | Only v3 is supported. The old v2 spelling is rejected, but the selector remains an unnecessary public surface to remove in a later input sweep. |
| hdf5 | Legacy optional backend. Shared checkpoints use native archives; HDF5 is outside the initial unification release. |
| subworlds/subworld_ranks and state_parallel* | Two existing scheduling interfaces. Shared native sessions now follow participating Worlds; both existing scheduler controls remain. |
| archive, fock_json_file, prefix | Older application-level path controls; workflow/reference ownership must be settled before deleting or relocating them. |
| seed.es_y, seed.freq_tol, dalton.dir | Active import behavior, distinct from the retired automatic stage-selection/warmup controls. |

Layout update, November 20: SCF/response single-point stores now live at the run
root, with separate analysis directories; see the
[current layout](../../src/apps/madqc/README.md#current-run-directory). Generations
still serve exact input bindings, atomic publication and retention of live
sources; attempts separate result source maps and iteration comparisons. Neither
is a scientific accuracy criterion. The layout migration did not change those
internal records or archive retention. A further simplification would need to
preserve dependency invalidation and atomic replacement.
