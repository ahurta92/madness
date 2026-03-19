# molresponse (Legacy) Notes

This directory contains the original `molresponse` implementation (`src/apps/molresponse/molresponse.cc` and related classes).

For the current application-interface workflow, prefer:

```bash
madqc --wf=response <input_file>
```

If a standalone executable is needed, use:

```bash
molresponse2 <input_file>
```

`molresponse2` is wired through the same response workflow setup as `madqc --wf=response` (via `src/madness/chem/WorkflowBuilders.hpp`), so response orchestration changes only need to be maintained in one place.

## Maintenance Guidance

1. Treat `madqc --wf=response` as the canonical entry point.
2. Keep `molresponse2` as a thin compatibility wrapper only.
3. Avoid adding new duplicated orchestration logic in `molresponse2.cpp` or `madqc.cpp` response branch.
