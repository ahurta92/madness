# Fresh named moldft → molresponse calculations

Run the checked-in H2 example with a new output directory:

```
MAD_NUM_THREADS=2 madqc named_h2.in --prefix=H2-run
```

`named_h2.in` defines a molecule, `moldft gs`, and `molresponse optical`.
The ground stage requests HF energy and dipole; the response stage requests
`zz` polarizability at frequencies `0.0` and `0.02` atomic units. These are
example numerical settings, not new solver defaults.

Named stages run in file order. Each `moldft` names an earlier molecule, and
each `molresponse` names an earlier `moldft` reference. Names are case-sensitive,
start with a letter, and otherwise contain letters, digits, `_`, or `-`.
All names must be unique. Each block ends with `end`; nested blocks are unsupported.
The complete deck is checked before any stage starts.

Every stage has its own typed parameters and directory. Two response stages
can reference the same ground stage while requesting different frequencies
and numerical settings. Geometry and the Hamiltonian come from that reference;
response `protocol`, `dconv`, and `maxiter` belong to the response stage. Omitted
response protocol/maxiter use the existing parameter-group defaults; omitted
response dconv uses the existing `100 * finest protocol threshold` derivation.
No response setting silently tightens the upstream calculation.

The first implementation supports these controls:

| Block | Controls |
|---|---|
| `molecule NAME` | Ordinary atom coordinates, `units`, `eprec`, `no_orient` |
| `moldft NAME` | `molecule`, `xc hf`, `protocol`, `dconv`, `econv`, `maxiter`, `L`, `k`, `charge`, `localize`, `print_level`, `maxsub`, `maxrotn`, `requested_properties [energy, dipole]` (either or both) |
| `molresponse NAME` | `reference`, `protocol`, `dconv`, `maxiter`, `kain`, `print_level`, `requested_properties [polarizability]`, `polarizability.frequencies`, `polarizability.frequency_units au`, `polarizability.components` |

Ground states must be closed shell. Response frequencies and components must
be explicit; components are Cartesian pairs such as `[zz]` or `[xx, xy]`.
The engine may calculate extra tensor elements needed to assemble the request;
`properties.polarizability.values` selects the requested components. The raw
engine evidence is retained alongside them. Response protocols finer than the
reference threshold/order are rejected. A response `k` override is unsupported.
Frequencies and response protocol rungs that collide at the engine's existing
filename-key precision are rejected. Unknown controls and unsupported properties
fail rather than being ignored.
Command-line scientific overrides and the old `--wf` option do not apply to
named decks; put settings in their owning block. Existing unnamed decks retain
their workflow interface.

Outputs go to `<prefix>/<stage>/`, without changing the process working
directory. `<prefix>/results.json` indexes the stage results; each stage also
has `results.json`. SCF's resolved input and engine files are under its stage,
and response metadata/state files are under the response stage. Results carry
method, units, effective reference cell, requested targets, and the existing
engine's convergence evidence. An unconverged stage is recorded and causes a
nonzero exit before descendants run.

Energy is in hartree, dipole in `e*a0`, and polarizability in `a0^3`; response
frequency is in atomic angular-frequency units. Property records identify their
fresh state and, for response, the source ground state. `fresh_run_state` IDs
are scoped by a unique run ID; `checkpoint_generation` is explicitly null.

This is milestone 2a: **fresh runs only**. An existing output root is refused.
The engine's native archives are intermediate output, not published generations
of the new shared checkpoint system. Ground-state checkpoint publication,
response bundle restart, reuse/tighter-target continuation, and the full measured
convergence contract are milestones 2b–2d. Beta, Raman, 2PA, excited states,
Nemo, MP2, CC2, and optimization are outside this slice.

Tests: `test_named_stages` covers block parsing, ownership and errors in serial
and two ranks. `test_named_stages_e2e.py` compares fresh H2 energy/alpha with the
existing workflow, exercises two response instances, preserves an existing output
root on refusal, and checks the cell after automatic centering of translated H2.
Its numerical runs retain logs and JSON beneath a temporary test directory.
On OpenMPI, assign sufficient CPUs per rank (for example `--bind-to none` for
local testing); do not run threaded tests with one core assigned per rank.
