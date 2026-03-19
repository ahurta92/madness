# Phase N Validation - <Title>

Date:

- <YYYY-MM-DD>

Related roadmap phase:

- <phase number and name>

Related notes:

- `<path to phase notes>`

## Validation Scope

- <what this validation is meant to prove>
- <what is explicitly not covered>

## Environment

- Repository root: `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next`
- Build directory: `<path>`
- Execution mode: `<serial | subgroup | mpi | interactive job>`
- Loader / environment script: `<load_40core.sh | load_xeonmax.sh | other>`
- Slurm job id if used: `<job id>`

## Commands Run

```bash
<command>
```

```bash
<command>
```

## Tests / Checks

| Check | Result | Notes |
| --- | --- | --- |
| <test> | <pass/fail> | <notes> |
| <test> | <pass/fail> | <notes> |

## Variant Coverage

| Variant | Covered in this validation? | Notes |
| --- | --- | --- |
| `static_restricted` | <yes/no> | <notes> |
| `static_unrestricted` | <yes/no> | <notes> |
| `dynamic_restricted` | <yes/no> | <notes> |
| `dynamic_unrestricted` | <yes/no> | <notes> |

## Metadata / Restart Checks

- `response_metadata.json` parity checked? `<yes/no>`
- Workflow metadata parity checked? `<yes/no>`
- Restart source verified? `<yes/no>`
- Restart support mode verified? `<yes/no>`

## Results Summary

- <result>
- <result>

## Remaining Validation Gaps

- <gap>
- <gap>
