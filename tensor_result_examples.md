# Grouped alpha and beta results

**Grouping and shared state-metric placement approved — 2026-11-11.** This
extends the earlier [output example](output_record_examples.md) with compact
state measurements in `results.json`. The 2026-11-12 consolidation aligns this
excerpt with the [record contract](design_draft.md#requested-properties-and-shared-state-metrics)
for review. Numbers remain invented. Version, attempt, `metric_definitions` and
`metric_set` are omitted from this excerpt; they are required in the full record.
The calculation's source-map-owning `attempt_id` inherits the enclosing results
run attempt only when equal; otherwise it is explicit on this calculation.
Individual numerical producing attempts remain in the source map.

Use static alpha and beta to keep the example small. Alpha has nine components;
beta has twenty-seven. Both use the same three dipole response states here, so
show their measurements once at the response-calculation level. Each property's
frequency or frequency combination identifies its scientific context; no list of
dependencies is repeated for every tensor component.

```json
{
  "calculation_id": "response",
  "application": "molresponse",
  "model": "TDHF",
  "frequency_units": "atomic_angular_frequency",
  "protocol": {"threshold": 1e-06, "k": 8},
  "properties": {
    "polarizability": {
      "units": "bohr^3",
      "values": [
        {
          "frequency": 0.0,
          "status": "complete",
          "components": {"xx": 10.0, "xy": 0.2, "xz": 0.3, "yx": 0.2, "yy": 11.0, "yz": 0.4, "zx": 0.3, "zy": 0.4, "zz": 12.0},
          "availability": {}
        }
      ],
      "value_type": "real"
    },
    "hyperpolarizability": {
      "units": "e^3*bohr^3/hartree^2",
      "values": [
        {
          "process": "static",
          "frequencies": {"B": 0.0, "C": 0.0},
          "output_frequency": 0.0,
          "status": "complete",
          "components": {
            "xxx": 1.1,
            "xxy": 0.2,
            "xxz": 0.3,
            "xyx": 0.2,
            "xyy": 0.4,
            "xyz": 0.5,
            "xzx": 0.3,
            "xzy": 0.5,
            "xzz": 0.6,
            "yxx": 0.2,
            "yxy": 0.4,
            "yxz": 0.5,
            "yyx": 0.4,
            "yyy": 2.2,
            "yyz": 0.7,
            "yzx": 0.5,
            "yzy": 0.7,
            "yzz": 0.8,
            "zxx": 0.3,
            "zxy": 0.5,
            "zxz": 0.6,
            "zyx": 0.5,
            "zyy": 0.7,
            "zyz": 0.8,
            "zzx": 0.6,
            "zzy": 0.8,
            "zzz": 3.3
          },
          "availability": {}
        }
      ],
      "value_type": "real"
    }
  },
  "state_metrics": [
    {
      "frequency": 0.0,
      "states": {
        "x": {
          "values": {"bsh_residual": 3e-07, "density_change": 4e-07, "kinetic_change": 6e-08},
          "availability": {"energy_change": {"status": "not_applicable", "reason": "fixed_driving_frequency"}}
        },
        "y": {
          "values": {"bsh_residual": 4e-07, "density_change": 5e-07, "kinetic_change": 7e-08},
          "availability": {"energy_change": {"status": "not_applicable", "reason": "fixed_driving_frequency"}}
        },
        "z": {
          "values": {"bsh_residual": 5e-07, "density_change": 6e-07, "kinetic_change": 8e-08},
          "availability": {"energy_change": {"status": "not_applicable", "reason": "fixed_driving_frequency"}}
        }
      }
    }
  ]
}
```

The metrics are a snapshot of the exact states used for these properties.
Ground-state measurements appear with the SCF calculation. Shared report context
supplies metric definitions and units; the full checkpoint retains comparison
identities, relative scales and publication/source identities. This small
results excerpt omits those surrounding records. No driven energy-change metric
is manufactured, and direct VBC tasks do not acquire residuals.

A full tensor request accounts for all nine or twenty-seven component names.
A partial request has fewer names. Unavailable requested components move to the
availability map with a reason; unrequested components stay absent. The equal
values above are illustrative, not an instruction to fill entries by symmetry.
Declared complex quantities retain the agreed `[real, imaginary]` encoding.

For multiple frequencies, add records under each property's `values` and under
`state_metrics`. Beta retains its process and ordered B/C frequencies; its
measurements must cover the actual contributing states at every required
frequency. A single x/y/z group works here because the example is static. For
other operators or distinct states at the same frequency, use unambiguous
application-assigned state labels; x/y/z alone is not a universal identity.

The approved layout reports shared state measurements once per response calculation,
which avoids repeating them in both alpha and beta. Published results bind one
consistent set of source states; scientific reuse and compatibility still use
the exact runtime/checkpoint evidence. Metrics describe input-state convergence,
not a numerical error bound on the tensor.

The beta process and dimensional convention follow the existing
[response formatter](src/apps/madqc/ResponseProperties.hpp) and the maintained
[static beta deck](src/examples/qc/response_lih_beta/response_lih_beta.in).
This example does not add unsupported response processes or change the agreed
checkpoint, timing, or convergence policy.
