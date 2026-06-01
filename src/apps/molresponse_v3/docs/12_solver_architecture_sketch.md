# Solver architecture sketch (v3)

Status: **sketch**, not yet wired. Captures the design that emerged from
the 2026-05-15 TDA debugging session in `src/apps/molresponse/` — see
the closing notes in this file for the specific bugs that motivate this
architecture.

## Problem statement

The unified solver in `src/apps/molresponse/` tries to handle TDA, RPA,
static, frequency-dependent, closed-shell, and (eventually) open-shell
all in the same call paths. Each operation does a runtime branch on
`calc_type == "tda" / "full" / "static"` or `r_params.tda()`, and a
single `X_space` carries both `.x` and `.y` regardless of whether
`.y` is meaningful for the calculation.

Two classes of bug followed from that unification:

1. **Empty-block UB.** Closed-shell TDA leaves `.y` blocks empty and
   never populated; the bare `response_space(world, m, n)` ctor leaves
   inner vectors empty by default. `X_space::to_vector` then indexes
   `y[ai][j]` out-of-bounds during ordinary `+` / `-` / `*` ops,
   silently corrupting results. This produced a factor-of-2 error in
   `<X|Gamma|X>` for closed-shell TDA that was invisible during RPA
   testing (where `.y` is populated).
2. **Aliased state.** `ResponseBase::Chi` is a class member but also
   appears as a function parameter `chi` in many subroutines. Multiple
   sites read `Chi` (member, possibly pre-rotation) when they should
   have used the parameter (post-rotation). Found via grep, not by
   tests.

## Design principles

1. **One solver class per (response type × shell type).** No runtime
   branching on `calc_type`. Compile-time errors instead of silent
   data-dependent UB.

2. **State is explicit and immutable across step boundaries.** Each
   solver step takes `IterationState in` and returns `IterationState
   out`. No class-member trial vectors. The data dependencies of every
   subroutine are visible in its signature.

3. **Sequence control is a free template, not a method.** A single
   `iterate(solver, state, policy)` function drives every solver. The
   exact step order is readable in one place; specific solvers
   implement the steps but do not choose the order.

4. **No god-class.** `ResponseBase` ownership of `Chi`, `rho_omega`,
   `omega`, etc. goes away. The solver class owns *policy* (KAIN
   subspace size, convergence thresholds, step-restriction max) and
   reusable handles (ground orbitals, BSH operators); iteration data
   flows through `IterationState`.

5. **No data-type-mixing inside a single function.** A function that
   takes a closed-shell response vector cannot accidentally pull in
   `.y` because the type doesn't have one.

## Types

### `ResponseTarget<Shell, Mode>`

Invariant under iteration. Describes *what response we are computing*.

```cpp
struct ClosedShell {};
struct OpenShell  {};

struct ExcitedState {};                       // many roots, omegas to solve for
struct FrequencyDependent { double omega; };  // one omega per perturbation, given

template <typename Shell, typename Mode>
struct ResponseTarget;

template <> struct ResponseTarget<ClosedShell, ExcitedState> {
  std::vector<real_function_3d> phi0;     // ground occupied orbitals
  Tensor<double> ground_energies;
  int n_roots;                            // how many excited states
  // ... XC functional handle, basis info, etc.
};

template <> struct ResponseTarget<OpenShell, ExcitedState> {
  std::vector<real_function_3d> phi0_alpha;
  std::vector<real_function_3d> phi0_beta;
  Tensor<double> ground_energies_alpha;
  Tensor<double> ground_energies_beta;
  int n_roots;
  // ...
};

template <> struct ResponseTarget<ClosedShell, FrequencyDependent> {
  std::vector<real_function_3d> phi0;
  Tensor<double> ground_energies;
  double omega;
  vector_real_function_3d perturbation;   // dipole, e.g.
  // ...
};
```

Restart artifact: serialize `ResponseTarget` once at job start, never
again.

### `IterationState<Shell, Mode>`

Volatile. Carries the *current iterate and where the solver is in its
convergence trajectory*.

```cpp
template <typename Shell, typename Mode>
struct IterationState;

template <> struct IterationState<ClosedShell, ExcitedState> {
  response_space x;          // size [n_roots × n_occ], no Y block at all
  Tensor<double> omega;      // current eigenvalues
  // diagnostics / convergence
  X_space last_residual;
  double last_bsh_residual;
  double last_density_residual;
  int    iter;
  double protocol_thresh;
  // restart cursor
  std::optional<KainHistory> kain_history;
};

template <> struct IterationState<ClosedShell, FrequencyDependent> {
  response_space x;          // size [n_perturbations × n_occ]
  response_space y;          // size [n_perturbations × n_occ] — real, not zero
  // ... rest analogous
};
```

Key point: for closed-shell TDA, `IterationState` has **no `.y` field
at all**. The type system prevents anyone from writing
`state.y[b][p] = ...` and prevents `to_vector` from emitting a `2*n*m`
flat vector. The UB we just fixed becomes uncompilable.

Restart artifact: snapshot `IterationState` at each protocol boundary.

## Solver classes

One concrete solver per `(Shell, Mode)`. They are duck-typed: each
exposes the same method names so a common template can drive them, but
they are not in an inheritance hierarchy.

```cpp
class ClosedShellTdaSolver {
public:
  using State  = IterationState<ClosedShell, ExcitedState>;
  using Target = ResponseTarget<ClosedShell,  ExcitedState>;

  ClosedShellTdaSolver(const Target& target, SolverPolicy policy);

  // The six steps every TDA / RPA / FD solver implements. Each is a
  // pure (or close to pure) transformation of the iteration state.
  // No class-member trial data is read or written.
  State compute_lambda(const State& in) const;
  State rotate          (const State& in) const;
  State compute_theta   (const State& in) const;
  State bsh_update      (const State& in) const;
  State compute_residual(const State& in) const;
  State kain_step       (State in);       // owns its KAIN history; needs &
  State step_restrict   (const State& in) const;

  // Policy queries — drives whether iterate() calls a given step.
  bool use_kain() const;
  bool use_step_restriction() const;

  // Convergence + serialization
  bool   converged(const State& s) const;
  void   save     (const State& s, const std::string& path) const;
  State  load     (const std::string& path) const;

private:
  const Target& target_;
  SolverPolicy policy_;
  KainHistory kain_;                          // owned by solver, not state
};

class ClosedShellRpaSolver {
  using State  = IterationState<ClosedShell, ExcitedState>;
  // Same method names, but State has .y and the bodies handle it.
};

class OpenShellTdaSolver {
  using State  = IterationState<OpenShell, ExcitedState>;
  // Same method names, different math under the hood.
};

class ClosedShellFrequencyDependentSolver { /* ... */ };
// ... etc.
```

## The driver

The iteration loop is **one function, written once**. It takes any
solver and any compatible state and runs the standard sequence. No
solver controls the order; no order is hidden inside a solver class.

```cpp
template <typename Solver>
auto iterate(Solver& solver, typename Solver::State state,
             const ConvergencePolicy& cp) -> typename Solver::State {
  for (int k = 0; k < cp.max_iters; ++k) {
    state = solver.compute_lambda(state);
    state = solver.rotate(state);
    state = solver.compute_theta(state);
    state = solver.bsh_update(state);
    state = solver.compute_residual(state);

    if (solver.use_kain()) {
      state = solver.kain_step(state);
    }
    if (solver.use_step_restriction()) {
      state = solver.step_restrict(state);
    }

    if (solver.converged(state)) break;
  }
  return state;
}
```

User-level call site:

```cpp
auto target = build_closed_shell_tda_target(world, ground_state, r_params);
auto state0 = make_initial_guess(target, restart_path);
ClosedShellTdaSolver solver(target, policy);
auto state_final = iterate(solver, state0, cp);
solver.save(state_final, save_path);
```

## What this prevents

- **`.y` UB for closed-shell**: closed-shell `IterationState` has no
  `.y`, so it cannot be accessed empty. The `response_space::operator+`
  bug we hit on 2026-05-15 (empty inner vectors → OOB in `to_vector` →
  factor-of-2 in `<X|Gamma|X>`) becomes a compile error in v3.

- **Aliased trial state**: there is no `Solver::Chi_` member to be
  stale or out of sync with a parameter. Every method takes `State`
  by const-ref and returns a new `State`. Data flow is explicit.

- **Cross-mode contamination**: a `compute_lambda(State<ClosedShell, ExcitedState>)`
  cannot accidentally call into a frequency-dependent specialization;
  the type doesn't match.

- **Hidden step-order changes**: the order of operations is one
  template function. Reviewing a change to it is a single-file diff,
  not a hunt across every solver.

## What this still has to figure out

1. **Shared helpers** that don't care about shell or mode (BSH operator
   build, projection, plain dot/inner). Top-level free functions or a
   shared `ResponseKernel` utility — TBD. Either way, they take their
   inputs by ref and don't read class members.

2. **KAIN state lifetime**. KAIN history naturally lives in the solver
   class (it's policy + history, not iteration data). `kain_step`
   needs the solver as `&`, not `const&`. The driver template handles
   this via `Solver& solver`.

3. **Macrotask integration**. The per-state loops inside
   `compute_lambda` / `compute_theta` etc. should be amenable to
   replacement with a MacroTask dispatch when distributed perf
   demands it. The signature `State f(const State&) const` doesn't
   constrain how the body is implemented.

4. **Sharing math across closed-shell-TDA and closed-shell-RPA**.
   Bodies differ but some core operators don't. Don't reach for
   inheritance — when you find yourself duplicating, extract a free
   function or a small struct. Inheritance hierarchies across these
   specializations have historically grown messy.

## Why this is worth the refactor

Compare current `src/apps/molresponse/` `compute_gamma_tda`:

- 200+ lines
- Reads class members `ground_orbitals`, `shared_coulomb_operator`
- Mutates state via `gamma.x = (J*2) - k1_x*c + W` where each operand
  has subtle invariants that determine whether the result is correct
- Internal `orbital_load_balance` mutates the global function pmap

vs. what it would look like in v3:

```cpp
IterationState<ClosedShell, ExcitedState>
ClosedShellTdaSolver::compute_theta(const State& in) const {
  State out = in;
  out.theta = make_theta_closed_shell_tda(
      target_.phi0, in.x, target_.xc_operator);
  return out;
}
```

The mathematical content lives in `make_theta_closed_shell_tda`, which
takes its inputs by ref, returns a value, and has no shell- or
mode-branching inside. Auditable in 30 lines.

## Bugs this session that motivated this design

1. **`response_space::operator+` constructed its result with the bare
   ctor**, which left inner vectors empty. `from_vector` then wrote
   to OOB indices and downstream arithmetic mixed in garbage. For
   closed-shell TDA this manifested as `<X|Gamma|X>` being roughly
   half of `2J - K1`, which produced the wrong Lambda eigenvalues
   and a wrong converged answer.

2. **`X_space` constructors leave `.y` empty** for closed-shell TDA;
   `X_space::to_vector` / `from_vector` pack `(X, Y)` interleaved and
   read OOB on `.y`. Every `X_space` arithmetic op was subtly wrong.

3. **`Chi.x` vs `chi.x`** in `bsh_update_excited`: function used the
   class member instead of the parameter; for the bound-state regime
   the shift terms were zero so the bug was silent. Fixed by adding
   an explicit `chi` parameter.

4. **Lambda / theta truncation discipline** in `compute_response_potentials`:
   legacy truncates after assembly, current didn't. With v3's
   "every step returns a new State", the truncation can be a single
   place at the end of each step.

5. **Unified `rotate_excited_space`** routed TDA through indefinite-metric
   sygvp + cluster polar-decomp that was designed for RPA's eigenvalue
   structure; for TDA's tightly-clustered bound-state omegas this
   could scramble state identity iter-to-iter. Fixed by splitting
   the function into explicit TDA / RPA branches.

All five of these are class-of-bug eliminated by the v3 typing.
