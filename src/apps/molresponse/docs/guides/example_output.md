# Example output

What a response calculation actually prints and writes. Every excerpt below is
**real output from a run of this code**, lightly trimmed for length — not
illustrative mock-ups.

Two output surfaces:

- **`<prefix>.out`** and **`<run>/summary.txt`** — the human-readable summary.
- **`<run>/results.json`** — the machine-readable record: every property, each
  item identified by its fields and carrying its status and convergence. This is
  what you script against; see [results_file.md](results_file.md).

## A polarizability run, end to end

Deck (water, HF, static α):

```text
dft
    xc        hf
    protocol  [1e-4]
    l         200.0
    prefix    "mad"
end
response
    dipole             true
    dipole.directions  xyz
    dipole.frequencies [0.0]
end
molecule
    units atomic
    no_orient true
    O    0.00000000   -0.22490589    0.00000000
    H    1.45235000    0.89962300    0.00000000
    H   -1.45235000    0.89962300    0.00000000
end
```

`madqc --wf=response h2o.in` → `h2o.out`:

```text
======================================================================
  MADQC RESULTS SUMMARY
======================================================================

  Task 0 : SCF  (model = scf)
  ----------------------------------------------------------------------
    Molecule         : H2O
    Total energy     :      -76.065377246 Ha       -2069.8444 eV
    Orbital eps (a)  : -20.5621 -0.9593 -0.9593 -0.7225 -0.5092  (Ha)
    Converged        :  thresh = 0.0001  dconv = 0.0001
    Wall time        :  0.0 s   (1 MPI x 1 threads)

  Task 1 : RESPONSE
  ----------------------------------------------------------------------
    alpha[1e-04_k6](w=0.000)
      x :      9.482637     -0.000000     -0.000000
      y :     -0.000000      8.714917      0.000000
      z :      0.000000      0.000000      7.973547

======================================================================
```

Note the shape: the workflow ran the ground state as task 0 and the response as
task 1, and the α label carries the **protocol key** (`1e-04_k6`) and the
frequency, because both are part of the result's identity. The off-diagonal
elements are zero to ~10⁻⁷ here, which is the symmetry of this geometry showing
up as a numerical check for free.

## Excited states

Per converged root, the solver reports the excitation energy and its transition
properties (real output, water, root 0):

```text
   Response Function 0		8.54365210 eV
   --------------------------------------------

   Transition Dipole Moments
   X: 0.00000000   Y: 0.00000000   Z: -0.45734020

   Dipole Oscillator Strength: 0.04378044

   Transition Quadrupole Moments
                  X                Y                Z
   X       0.00000000      -0.00000000      -0.00000187
   Y      -0.00000000       0.00000000       0.24705616
   Z      -0.00000187       0.24705616      -0.00000000

   Dominant Contributions:
                  x          y
   Occupied 4   0.99991358 0.02164660
   Occupied 3   0.02094647 0.02793969
   Occupied 1   0.01199088 0.00997191
```

**Dominant Contributions** is the useful diagnostic for identifying *which*
state you got: it lists the occupied orbitals carrying the excitation, with
separate `x` (excitation) and `y` (de-excitation) weights. Here the state is
essentially a single excitation out of occupied orbital 4 (weight 0.9999) with
negligible de-excitation character — as expected for a low valence excitation.
Use this, not the root index, when matching states across calculations.

## Two-photon absorption

For a 2PA run the engine prints the transition tensor per state and then the
rotationally-invariant observables (real output, LiH — note roots 2 and 3 are a
degenerate pair):

```text
                  | Two-photon transition tensor S |
                  +--------------------------------+
     -----------------------------------------------------------------
      No  Energy(eV)     Sxx     Syy     Szz     Sxy     Sxz     Syz
     -----------------------------------------------------------------
       1       4.06    25.792  25.792  75.255   0.000  -0.000   0.000
       2       5.07    -0.000   0.000   0.000  -0.000  -7.372 -42.533
       3       5.07    -0.000  -0.000  -0.000   0.000 -42.533   7.372
     -----------------------------------------------------------------

     D = 2*Df+4*Dg (Linear);  D = -2*Df+6*Dg (Circular)
     Df = sum_ij S_ii*S_jj /30;   Dg = sum_ij S_ij^2 /30
     sigma = D*(E/2)^2*AU_TO_GM (GM, 0.1 eV FWHM);  R = (-Df+3Dg)/(Df+2Dg)
```

The definitions are printed with the numbers deliberately, so a table pasted into
a notebook is self-documenting. **Read the invariants, not the components,** for
degenerate states: roots 2 and 3 above distribute the same physics differently
between `Sxz` and `Syz` (the manifold is free to rotate), while their invariants
are identical — which is why any cross-code or cross-basis comparison must be
made on `Df`/`Dg`/`D`/`R`.

## The machine-readable record

`<run>/results.json` holds the same results, one entry per application. From
the `response_h2_es_rpa` example (H₂, HF, protocol 1e-4), a polarizability item
of the `response` entry:

```json
"polarizability": {
  "units": "bohr^3",
  "values": [
    { "frequency": 0.0, "status": "complete",
      "components": { "zz": 6.436397286910804 },
      "convergence": { "bsh_residual": 2.25e-06, "density_change": 3.42e-05 } }
  ]
}
```

and the state it was computed from, in the entry's `state_metrics`:

```json
"dipole_z-f0.00000": {
  "protocol": { "k": 6, "threshold": 0.0001 },
  "updates": 6,
  "status": "complete",
  "subjects": { "response": { "values": { "bsh_residual": 2.250088916218293e-06,
                                          "density_change": 3.4246290835799824e-05,
                                          "kinetic_change": 8.82517065159405e-06 } } }
}
```

Three conventions worth knowing:

- An entry holds the final level of the request; each state's earlier levels
  are in its node's `history.jsonl`.
- Items are identified by their fields (`frequency`; `frequency_b` and
  `frequency_c`; `root`), not by position.
- Every value carries `status` and a `convergence`, so a consumer can tell how
  trustworthy a number is without re-reading the log.

## Where files land

```
<prefix>/                          the run folder
  results.json                     the record above
  summary.txt                      results.json for people
  timings.json                     where the time went
  states/ground_state/             the ground state (archive, description.json, history.jsonl)
  states/response/<node>/          each response state, e.g. dipole_z-f0.00000, excited-full
  analysis/response/es_analysis__<protocol>.json   transition analysis of the excited states
<prefix>.out                       the human summary
<prefix>.calc_info.json            the workflow-level record (all tasks)
```

Response states are saved as node folders and are **restartable**: a re-run
reuses a state whose stored metrics meet the request, and otherwise resumes or
climbs from what is saved rather than starting over. See
[`madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md) for the
full deck reference and [`madqc/HDF5_IO.md`](../../../madqc/HDF5_IO.md) if you
want these archives as HDF5 instead of the native format.
