# molresponse_v2 Execution Notes

## Canonical Way To Run Response

Use `madqc` with the response workflow:

```bash
madqc --wf=response [options] [input_file]
```

This is the primary interface and should be used for new workflows, tests, and docs.

## Status Of `molresponse2`

`molresponse2` is now a compatibility wrapper.
It no longer owns solver workflow logic.
It dispatches to the same workflow construction used by `madqc --wf=response`.

```bash
molresponse2 [options] [input_file]
```

`molresponse2` only supports `response` workflow behavior.  
For other workflows (`scf`, `nemo`, `cc2`, `cis`, `oep`) use `madqc`.

## Single Source Of Workflow Logic

Response workflow wiring for both entry points is centralized in:

- `src/madness/chem/WorkflowBuilders.hpp`

Response implementation is centralized in:

- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/molresponse_v2/*` (solver/property/state machinery)

## Maintenance Rule

When changing response execution flow:

1. Update `WorkflowBuilders.hpp` and/or `MolresponseLib.hpp`.
2. Do **not** add standalone solver flow back into `molresponse2.cpp`.
3. Validate parity by running:
   - `madqc --wf=response ...`
   - `molresponse2 ...`
   and comparing generated `*.calc_info.json`.
