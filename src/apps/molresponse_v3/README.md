# molresponse_v3

Ground-up reimplementation of the MADNESS molecular response property
pipeline, designed for clarity, testability, and clean separation of
concerns.

## Status

Increment 0: skeleton app. Loads a ground-state checkpoint and prints
basic information.

## Relationship to other response codes

- `molresponse_legacy`: frozen reference for legacy excited-state behavior
- `molresponse_v2`: current production response code
- `molresponse_v3`: this code — incremental rebuild guided by design docs in `docs/`

v3 links against MADresponse2 during development to reuse working pieces.
As v3 matures, those dependencies are replaced. v3 is validated against
v2 at every increment.

## Design Documents

See `docs/` for the full architectural design.

## Build

From the MADNESS build directory:

```bash
cmake --build . --target molresponse_v3
```

## Run

```bash
molresponse_v3 --archive=<path_to_restartdata>
```

Or with MPI:

```bash
mpirun -n 1 molresponse_v3 --archive=<path_to_restartdata>
```
