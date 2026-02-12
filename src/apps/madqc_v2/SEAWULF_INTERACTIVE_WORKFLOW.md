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
