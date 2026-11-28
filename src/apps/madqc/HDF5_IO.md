# HDF5 I/O (optional)

madqc and the molresponse engine can store their **response states and restart
checkpoints** as HDF5 instead of the native MADNESS archive. This is an **opt-in**
feature: the native archive is the default and is always available, and a build
without HDF5 is byte-for-byte identical to one with the feature compiled but
unused.

## When you'd want it

- Interoperability / inspection with standard HDF5 tooling (`h5ls`, `h5py`, …).
- A single self-describing checkpoint blob rather than the native per-part
  archive family.

For ordinary runs the native backend is fine — you do not need HDF5.

## Enabling it

**1. Build with HDF5 (default OFF):**

```
cmake -DMADNESS_ENABLE_HDF5=ON …    # requires libhdf5; guarded by MADNESS_HAS_HDF5
```

All HDF5 code is behind that guard, so a normal (`OFF`) build carries none of it.

**2. Turn it on at runtime**, either scope:

```text
response
    hdf5  true         # response states only
end
```

or run-wide for every app that supports it:

```text
io
    backend  hdf5      # default: native
    nio      1         # number of I/O aggregator ranks
end
```

A stray `response.hdf5` / `io backend hdf5` in a deck for an app that does **not**
support HDF5 (moldft, cc2, cis, oep) is ignored — it neither flips the global flag
nor mis-stamps the backend.

## How it behaves

- The backend a launch writes with is recorded in `results.json` as
  `provenance.archive_io` (`native` or `hdf5`).
- HDF5 writes are a rank-0 gather at `nio=1`, covering closed-shell response
  states and checkpoints.
- Native and HDF5 states can coexist in one calc directory across runs; each
  loads correctly from its recorded backend.

## Not enabling it

Do nothing. The native MADNESS archive is the default; every feature works
identically, and restart/checkpoint files are the native per-part archives.
