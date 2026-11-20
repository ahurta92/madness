# Response restart after SCF refinement

Status: author clarification, November 19, 2026. Design only; runtime unchanged.
This replaces the earlier proposal to make general orbital alignment and old-SCF
archive retention part of ordinary response restart.

## Ordinary restart

Response starts from the converged occupied SCF orbitals. Its functions are paired
with that occupied-orbital representation. Tightening convergence tolerances does
not by itself change that pairing or require alignment.

| Change | Response behavior |
|---|---|
| Tighter response request; occupied SCF functions unchanged | Keep the saved functions and compare valid saved measurements with the current request. Reuse immediately if they qualify; otherwise continue from them. |
| Occupied SCF functions refined in the same model, ordering and localization representation | Continue from the saved paired response functions at their requested saved protocol. Rebuild affected operators/sources and obtain fresh response measurements against the refined SCF functions. |
| Explicit orbital rotation or changed localization representation | The application handles the known representation change explicitly before using the old response pairing. This is outside the ordinary tighter-tolerance path. |
| Imported basis-set seed with a different occupied-orbital representation | Use the import path's existing transformation/alignment. |
| Incompatible model, geometry or cell | Report the mismatch; do not infer an automatic scientific conversion. |

A valid starting point is distinct from a solution that already meets the new
request. If the ground-state functions or relevant equations actually change,
old response convergence and property values do not certify the changed problem.
If only requested tolerances change, valid saved measurements remain useful and
can avoid another evaluation. A new archive identifier is not itself evidence
that occupied orbitals rotated.

Keep per-state protocol progress. Reset KAIN and unavailable comparison history
when iteration resumes, and rebuild operators for the selected protocol. A saved
protocol absent from the request follows the already-agreed supported-conversion
rule. These rules need no additional user knobs.

## Application responsibility

SCF continuation should preserve occupied-orbital ordering and the representation
used to pair response functions, or make an actual rotation explicit. Molresponse
owns any corresponding response transformation. The shared checkpoint/controller
layer must not infer orbital matching or automatically maintain old orbital
frames for every numerical refinement.

Ordinary restart does not introduce a requirement to retain an old SCF archive
for alignment. Keep the intended one-current-checkpoint direction and normal
safe replacement. Simplifying the existing generation/retention implementation
remains separate work; this note does not claim that code has already changed.

The existing DALTON alignment routine remains relevant to imported seeds, whose
orbital representation is not controlled by the continued SCF calculation. Its
import-specific overlap cutoff and diagnostics are not a new restart policy.

## What the code already provides

- `solvers/shared_response_checkpoint.hpp`: native response loading and threshold/order
  conversion. A changed same-model SCF source retains response functions and saved
  protocol progress as a guess, with fresh measurements. Ordinary model changes
  are rejected before replacing response provenance.
- `solvers/fd_solver.hpp`: fresh measurements against current ground-state
  operators, and reset of KAIN/comparison history.
- `src/madness/chem/SCF.cc::get_fock_transformation`: calls
  `Localizer::undo_reordering` and `undo_degenerate_rotations`; ordinary canonical
  iteration uses this helper through `diag_fock_matrix`.
- `src/madness/chem/SCF.cc::prepare_final_orbitals`: source inspection confirms
  the canonical endpoint bypasses those protections by calling `sygvp` directly.
  The localized endpoint transforms only virtuals and leaves occupied functions
  untouched. A restart failure caused by this difference has not been reproduced.
- `tools/dalton_mra.hpp::occupied_gauge_rotation` and `solvers/dalton_import.hpp`:
  explicit imported-orbital alignment followed by occupied-space projection.

Paths without `src/` above are relative to `src/apps/molresponse/`.

## Scope after walkthrough

The author set aside the proposed final-canonical-preparation change. Keep SCF's
existing preparation behavior; neither that change nor a new endpoint continuity
test is a prerequisite for ordinary response restart. The existing evaluation of
the functions actually saved remains part of the shared convergence contract.

## Next integration boundary

Apply the same simple model-versus-numerical distinction to ordinary SCF and
response restart:

- Reject a model mismatch in the same calculation, including changed geometry
  or cell. SCF's current `plan()` instead permits changed scientific settings as
  a guess when the archive layout is loadable; geometry therefore needs an
  explicit model check before that fallback.
- Allow same-model numerical refinement. Response's current constructor rejects
  any changed exact SCF source before it can plan continuation; replace that
  blanket rejection with the agreed model check and fresh response measurements
  when SCF functions or relevant numerical operators changed.
- Keep saved-metric qualification when only the requested tolerances change.
  Iteration limits, reporting and checkpoint cadence do not change the model.

The existing settings snapshots contain both model and numerical details, so
comparing the entire snapshot is not the model check. Keep recording both, but
classify their differences explicitly. Deliberate seeds for another calculation
and workflow-managed geometry changes remain separate from ordinary restart;
optimizer restart remains deferred. These are pending integration changes.

## Focused implementation checks

1. Tighter response criteria with unchanged SCF: reuse qualifying evidence and
   iterate only states that need more work; do not rotate paired functions.
2. Same-model SCF refinement: preserve orbital pairing and response protocol
   progress, then check the continued response against the updated SCF functions.
3. Keep localization changes and DALTON seed transformations explicit and separate
   from ordinary restart. Preserve their existing scientific regression coverage.
4. Reject real model mismatches and prevent stale response/property completion
   from being reported after the ground-state functions actually change.

No new calculation, performance claim or runtime change accompanies this note.
