# Personal Seawulf Interactive Workflow (MADQC)

This is a personal runbook for interactive development and testing on Seawulf HBM nodes.
It is intentionally separate from `AGENTS.md`.

## 1) Allocate an HBM interactive node

Use 1 node and 8 MPI tasks mapped to 8 NUMA regions.
Do not set `--cpus-per-task`.

```bash
salloc -p hbm-short-96core --nodes=1 --ntasks-per-node=8 --time=01:00:00
```

## 2) Load the build/run environment

Use the same environment as the CMake kit in VSCode.
From shell, source whichever script matches the build:

```bash
source /path/to/load_xeonmax.sh
# or
source /path/to/load_hdf5.sh
```

## 3) Build `madqc`

From repository root:

```bash
cmake --build build --target madqc -j 16
```

## 4) Run response test workload with NUMA mapping

Use the runtime pattern below for quick interactive checks:

```bash
MAD_NUM_THREADS=10 \
mpirun --map-by numa numactl --preferred-many=8-15 \
./build/src/apps/madqc_v2/madqc \
  --wf=response \
  --input=src/apps/madqc_v2/test_molresponse_h2o_alpha_beta_z.in
```

## 5) Regenerate molresponse reference JSON (when intentionally updating expected values)

Run with the test prefix used by the scripted test:

```bash
MAD_NUM_THREADS=10 \
mpirun --map-by numa numactl --preferred-many=8-15 \
./build/src/apps/madqc_v2/madqc \
  --wf=response \
  --prefix=mad_madqc_test_molresponse_h2o_alpha_beta_z.py \
  src/apps/madqc_v2/test_molresponse_h2o_alpha_beta_z.in
```

This produces:

```text
mad_madqc_test_molresponse_h2o_alpha_beta_z.py.calc_info.json
```

If the new output is correct, update the checked-in reference:

```bash
cp mad_madqc_test_molresponse_h2o_alpha_beta_z.py.calc_info.json \
  src/apps/madqc_v2/mad_madqc_test_molresponse_h2o_alpha_beta_z.py.calc_info.ref.json
```

## 6) Re-run the scripted regression

From `build/`:

```bash
ctest -R molresponse_h2o_alpha_beta_z --output-on-failure
```

For fast interactive development on Seawulf, keep the scripted test logic but
run it with your MPI/NUMA launcher via `MADQC_LAUNCHER`:

```bash
MADQC_LAUNCHER="mpirun --map-by numa numactl --preferred-many=8-15" \
MAD_NUM_THREADS=10 \
ctest -R molresponse_lih_alpha_raman_beta_xyz --output-on-failure
```

Default CI behavior is unchanged when `MADQC_LAUNCHER` is not set.

## 7) One-time geometry optimization for fixed-response tests

For response regression cases that should not optimize geometry during the test,
run a one-time SCF optimization in an interactive job and then copy the final
coordinates into the test input.

Important convergence note:
for `--optimize`, the geometry tolerances are derived from the **last**
`dft.protocol` value. Using `protocol=[1.e-4,1.e-6]` can make geometry
convergence too strict for quick test systems. Use a single-level protocol:
`protocol=[1.e-4]`.

Example:

```bash
MAD_NUM_THREADS=10 \
mpirun --map-by numa numactl --preferred-many=8-15 \
./build/src/apps/madqc_v2/madqc \
  --wf=scf \
  --optimize \
  --prefix=lih_opt_once \
  --dft="xc=hf; localize=new; maxiter=20; dconv=1.e-4; protocol=[1.e-4]; gmaxiter=20" \
  src/apps/madqc_v2/test_molresponse_lih_alpha_raman_beta_xyz.in
```

Then copy the final optimized coordinates from `lih_opt_once_opt.xyz` into the
test input geometry block so the scripted test remains a pure response run.
