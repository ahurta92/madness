#ifndef MOLRESPONSE_V3_CALC_CALC_MANAGER_HPP
#define MOLRESPONSE_V3_CALC_CALC_MANAGER_HPP

// ===========================================================================
// CalcManager — the scheduling core (doc 15, increment 15a).
//
// This header is the WORLD-FREE, PURE half of the calc manager: it turns a
// ResponsePlan (doc 14 Tier A) into a DAG of CalcNodes and emits, from the
// response session's plan for every state, the two-phase,
// perturbation-parallel schedule of (node, protocol) WorkItems.
//
// Nothing here touches madness::World, MPI, MRA functions, or the archive.
// Every function takes the aggregate metadata as a const nlohmann::json& so a
// test can pass a literal blob. The single World-bound seam is the abstract
// ICalcExecutor (declared at the bottom): the real FD/ES executor and
// CalcManager::run() live in calc/calc_executor.hpp (next increment) and are
// validated on the allocation.
//
// Design decisions vs doc 15:
//   * build_dag takes `int n_atoms`, not a full Molecule. The only reason the
//     planner needs the molecule is to expand the nuclear_all sentinel into
//     one node per atom; n_atoms is all that requires. CalcManager::build()
//     (executor side) extracts natom() from the Molecule and forwards it. This
//     keeps the whole header molecule-independent and trivially testable.
//   * planned_action / prerequisites_converged / schedule are keyed off the protocol step THRESHOLD
//     (a double), not a pre-built protocol_key string. k is derived once via
//     default_k_for_thresh, so callers only thread the ramp of thresholds.
//     (15a assumes the standard thresh->k table; an override_k ramp is a
//     later concern — see protocol_key_at().)
//   * prerequisites_converged resolves a node's prerequisites against the DAG (id -> node) so
//     it can reconstruct each prerequisite's metadata identity. In 15a no
//     concrete schedulable node carries a hard edge (DerivedFD "*" is a
//     symbolic expansion template, excluded from scheduling until ES
//     converges); prerequisites_converged becomes load-bearing for VBC in 15b and is
//     implemented + tested now so 15b inherits a working gate.
// ===========================================================================

#include "../Perturbations.hpp"          // Perturbation
#include "../ResponseProtocol.hpp"        // protocol_key, default_k_for_thresh
#include "../ResponsePropertyPlanner.hpp" // ResponsePlan, FDRequest, ESRequest, ...
#include "../solvers/response_metadata.hpp" // ResponseMetadata::freq_key

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Vocabulary (scheduler term -> response physics):
//   CalcNode       a response STATE to compute: x(ω) for a perturbation (FD),
//                  or an excited-state bundle (ES).
//   perturbation   the operator defining a unit of work (dipole_a, nuc_iα); the
//                  unit of parallelism — distinct perturbations are independent.
//                  (An ES bundle is its own unit, keyed by its id.)
//   seed_state     per perturbation, the state solved first (static ω=0 if
//                  requested, else lowest |ω|) that seeds frequency-continuation
//                  for the rest — solve α(0), continue to α(ω).
//   protocol step  one (thresh, k) MRA accuracy level on the coarse->fine ramp.
//   wave           a set of states with no prerequisite between them — safe to
//                  solve concurrently.
//   prerequisite   a state that must be CONVERGED before this one (correctness).
// ---------------------------------------------------------------------------

namespace molresponse_v3 {

// ---------------------------------------------------------------------------
// Schedulable unit.
// ---------------------------------------------------------------------------

enum class CalcKind {
  FD,         // concrete frequency-dependent dipole solve
  NuclearFD,  // nuclear-displacement FD (shares the fd_states subtree)
  ES,         // excited-state bundle
  DerivedFD,  // FD whose frequency comes from a converged ES root
  VBC,        // 15b: quadratic source (placeholder; not emitted yet)
};

/// One response state with its own protocol ladder. `id` is a stable string
/// that doubles as the dependency key (implicit-by-identity — see doc 15).
struct CalcNode {
  CalcKind            kind = CalcKind::FD;
  std::string         id;

  // FD / NuclearFD / DerivedFD payload
  Perturbation        pert;
  double              freq = 0.0;   // DerivedFD: 0 until expanded from a root

  // VBC payload: the C-side (B uses pert/freq above).
  Perturbation        pert_c;
  double              freq_c = 0.0;

  // ES payload
  bool                tda = true;
  int                 n_roots = 0;

  std::vector<double> protocols;    // coarse -> fine ladder (from the plan)

  // DerivedFD symbolic reference: "*" (all roots) or "es_root_0003".
  std::string         es_root_id;

  // DerivedFD: FD frequency = es_freq_factor * (ES root energy). 0.5 keeps the
  // FD off the response pole at ω = ωₙ (two-photon resonance, 2ω = ωₙ).
  // Default matches DerivedFDRequest (build_dag copies the request value in;
  // the two defaults are kept identical so a hand-built node can't silently
  // land on a different frequency rule). Unused on ordinary nodes.
  double              es_freq_factor = 0.5;

  // Hard edges (correctness): node ids that must be converged before this
  // node runs. Resolved by prerequisites_converged against the DAG + metadata.
  std::vector<std::string> prerequisites;

  /// A symbolic node is an expansion template, not directly schedulable:
  /// a DerivedFD that still references "*" (or no concrete frequency yet).
  bool is_symbolic() const {
    return kind == CalcKind::DerivedFD && es_root_id == "*";
  }
};

/// Reconcile verdict for one (node, protocol step). See the doc-15 table.
enum class NodeAction { Skip, Restart, Resume, Fresh };

/// Stringify a NodeAction for logs / the R1c scheduler trace (Output.diagnostics).
inline const char *node_action_name(NodeAction a) {
  switch (a) {
    case NodeAction::Skip:    return "skip";
    case NodeAction::Restart: return "restart";
    case NodeAction::Resume:  return "resume";
    case NodeAction::Fresh:   return "fresh";
  }
  return "?";
}

/// One (node, protocol) step of work.
struct WorkItem {
  const CalcNode *node = nullptr;
  double          thresh = 0.0;   // this protocol step's truncation threshold
  NodeAction      action = NodeAction::Fresh;
};

/// DESIGN B6, D-016: what the response session decided for one state. The
/// session gives every state a plan; schedule() refuses a state without one.
///   load    the saved state is used as it is: Skip at every rung. When it
///           meets today's request (`met`), a dependent state counts it as
///           converged.
///   resume  Skip every rung looser than start_thresh (a finished rung is never
///           redone). Each rung from start_thresh on is scheduled once in this
///           run; the session's node decides how far the state goes.
///   fresh   as resume, but the start rung is marked Fresh.
struct StatePlan {
  enum class Kind : std::uint8_t { load, resume, fresh };
  Kind        kind = Kind::fresh;
  double      start_thresh = 0.0;
  std::string loaded_key;   ///< load: protocol key of the saved state loaded
  bool        met = false;  ///< load: the saved state meets today's request
};
/// The plans of one run, by fd_node_id.
using StatePlans = std::map<std::string, StatePlan>;
/// Rungs this run has solved, as "<node id>|<protocol key>".
using SolvedRungs = std::set<std::string>;

// ---------------------------------------------------------------------------
// Identity helpers (free functions — also used by the executor later).
// ---------------------------------------------------------------------------

/// Canonical protocol_key for a protocol step, derived from the threshold alone using
/// the standard thresh->k table. (15a: override_k ramps unsupported.)
inline std::string protocol_key_at(double thresh) {
  return protocol_key(thresh, default_k_for_thresh(thresh));
}

inline std::string fd_node_id(const Perturbation &pert, double freq) {
  return "fd:" + pert.description() + "@" + ResponseMetadata::freq_key(freq);
}
inline std::string es_node_id(bool tda, int n_roots) {
  return std::string("es:") + (tda ? "tda" : "full") + "_n"
       + std::to_string(n_roots);
}
inline std::string derived_fd_node_id(const Perturbation &pert,
                                       const std::string &es_root_id) {
  return "dfd:" + pert.description() + "@" + es_root_id;
}
inline std::string vbc_node_id(const Perturbation &pb, const Perturbation &pc,
                               double fb, double fc) {
  return "vbc:" + pb.description() + "__" + pc.description() + "@"
       + ResponseMetadata::freq_key(fb) + "_" + ResponseMetadata::freq_key(fc);
}
/// The VBC source's node folder name: the parts of vbc_node_id joined by
/// hyphens. A description never holds a hyphen and each frequency part starts
/// with 'f', so two different sources never share a name.
inline std::string vbc_node_name(const Perturbation &pb, const Perturbation &pc,
                                 double fb, double fc) {
  return "vbc-" + pb.description() + "-" + pc.description() + "-"
       + ResponseMetadata::freq_key(fb) + "-" + ResponseMetadata::freq_key(fc);
}

/// The perturbation a node belongs to — the unit of parallelism (distinct
/// perturbations are independent and share every wave). FD / NuclearFD /
/// DerivedFD key on their perturbation operator; an ES bundle is its own unit,
/// so it keys on its own id.
inline std::string perturbation_key(const CalcNode &n) {
  if (n.kind == CalcKind::ES || n.kind == CalcKind::VBC) return n.id;
  return n.pert.description();
}

// ---------------------------------------------------------------------------
namespace detail_calc {

/// Sorted (coarse -> fine) deduplicated protocol ramp.
inline std::vector<double> sorted_ramp(std::vector<double> r) {
  std::sort(r.begin(), r.end(), std::greater<double>{});
  r.erase(std::unique(r.begin(), r.end(),
                      [](double a, double b) {
                        return std::abs(a - b) <= 1e-15 * std::max(a, b);
                      }),
          r.end());
  return r;
}

inline bool ramp_contains(const std::vector<double> &ramp, double thresh) {
  for (double t : ramp)
    if (std::abs(t - thresh) <= 1e-12 * std::max(t, thresh)) return true;
  return false;
}

inline bool entry_converged(const nlohmann::json &e) {
  return e.value("converged", false);
}

/// FD status entry at an exact (pert, key, freq) or nullptr. D-115: the
/// frequency key is "f%.5f", so two frequencies closer than its last digit
/// share it; an entry whose stored freq differs from `freq` beyond roundoff
/// holds another state and is not this one (it can still seed it, through
/// best_usable_fd_source_key). An entry with no stored freq (an older writer)
/// matches by its key, as before. A NaN `freq` skips the check.
inline const nlohmann::json *
fd_entry(const nlohmann::json &meta, const std::string &pert,
         const std::string &key, const std::string &fkey,
         double freq = std::numeric_limits<double>::quiet_NaN()) {
  if (!meta.contains("fd_states")) return nullptr;
  const auto &fs = meta["fd_states"];
  if (!fs.contains(pert) || !fs[pert].contains(key) ||
      !fs[pert][key].contains(fkey))
    return nullptr;
  const auto &e = fs[pert][key][fkey];
  if (freq == freq && e.contains("freq") && e["freq"].is_number() &&
      std::abs(e["freq"].get<double>() - freq) >
          1e-12 * std::max(1.0, std::abs(freq)))
    return nullptr;
  return &e;
}

/// Is this excited_states entry the bundle `node` asks for? The entry must
/// record the same type (tda/full) and root count. A legacy entry that records
/// neither is accepted (it predates the identity fields). Review finding C5: an
/// entry keyed on protocol alone let a converged TDA bundle satisfy a later RPA
/// or larger-n_roots request in the same calc dir.
inline bool es_entry_matches(const nlohmann::json &e, const CalcNode &node) {
  if (e.contains("type") &&
      e.value("type", std::string{}) != (node.tda ? "tda" : "full"))
    return false;
  if (e.contains("n_roots") && e.value("n_roots", 0) != node.n_roots)
    return false;
  return true;
}

/// ES bundle entry at a key that IS the requested bundle, or nullptr.
inline const nlohmann::json *es_entry(const nlohmann::json &meta,
                                      const std::string &key,
                                      const CalcNode &node) {
  if (!meta.contains("excited_states") ||
      !meta["excited_states"].contains(key))
    return nullptr;
  const auto &e = meta["excited_states"][key];
  return es_entry_matches(e, node) ? &e : nullptr;
}

} // namespace detail_calc

// ---------------------------------------------------------------------------
// build_dag — ResponsePlan -> nodes. Pure.
// ---------------------------------------------------------------------------

/// Expand a ResponsePlan into schedulable nodes:
///   * fd          -> one FD node each.
///   * es          -> one ES node each.
///   * derived_fd  -> one SYMBOLIC DerivedFD node each (es_root_id kept),
///                    prerequisites = {the plan's ES node id}. Expanded to
///                    concrete per-root nodes by the executor once the ES
///                    bundle converges (run(), next increment).
///   * nuclear_fd  -> all-atoms sentinel (pert.atom < 0) expands to one node
///                    per atom in [0, n_atoms); a concrete atom passes through.
///   * vbc         -> one node each; a C side with the all-atoms sentinel
///                    expands per atom the same way (full Raman).
std::vector<CalcNode> build_dag(const ResponsePlan &plan, int n_atoms) {
  std::vector<CalcNode> dag;

  for (const auto &r : plan.fd) {
    CalcNode n;
    n.kind      = CalcKind::FD;
    n.pert      = r.pert;
    n.freq      = r.freq;
    n.protocols = r.protocols;
    n.id        = fd_node_id(r.pert, r.freq);
    dag.push_back(std::move(n));
  }

  std::string first_es_id;
  for (const auto &r : plan.es) {
    CalcNode n;
    n.kind      = CalcKind::ES;
    n.tda       = r.tda;
    n.n_roots   = r.n_roots;
    n.protocols = r.protocols;
    n.id        = es_node_id(r.tda, r.n_roots);
    if (first_es_id.empty()) first_es_id = n.id;
    dag.push_back(std::move(n));
  }

  for (const auto &r : plan.derived_fd) {
    CalcNode n;
    n.kind       = CalcKind::DerivedFD;
    n.pert       = r.pert;
    n.es_root_id     = r.es_root_id;
    n.es_freq_factor = r.es_freq_factor;
    n.protocols      = r.protocols;
    n.id             = derived_fd_node_id(r.pert, r.es_root_id);
    if (!first_es_id.empty()) n.prerequisites.push_back(first_es_id);
    dag.push_back(std::move(n));
  }

  for (const auto &r : plan.nuclear_fd) {
    if (r.pert.is_all_atoms()) {
      for (int a = 0; a < n_atoms; ++a) {
        CalcNode n;
        n.kind      = CalcKind::NuclearFD;
        n.pert      = Perturbation::nuclear(a, r.pert.axis);
        n.freq      = r.freq;
        n.protocols = r.protocols;
        n.id        = fd_node_id(n.pert, r.freq);
        dag.push_back(std::move(n));
      }
    } else {
      CalcNode n;
      n.kind      = CalcKind::NuclearFD;
      n.pert      = r.pert;
      n.freq      = r.freq;
      n.protocols = r.protocols;
      n.id        = fd_node_id(r.pert, r.freq);
      dag.push_back(std::move(n));
    }
  }

  // One VBC node per (B,C) pair. A C side carrying the all-atoms nuclear
  // sentinel (full Raman) expands to one node per atom, exactly like the
  // nuclear_fd sentinel above, so its prerequisites name the concrete legs.
  for (const auto &r : expand_vbc_requests(plan.vbc, n_atoms)) {
    CalcNode n;
    n.kind          = CalcKind::VBC;
    n.pert          = r.pert_b;     // B side
    n.freq          = r.freq_b;
    n.pert_c        = r.pert_c;     // C side
    n.freq_c        = r.freq_c;
    n.protocols     = r.protocols;
    n.id            = vbc_node_id(r.pert_b, r.pert_c, r.freq_b, r.freq_c);
    // Hard prerequisites: the two converged first-order states (same protocol).
    n.prerequisites = {fd_node_id(r.pert_b, r.freq_b),
                       fd_node_id(r.pert_c, r.freq_c)};
    dag.push_back(std::move(n));
  }

  return dag;
}

// ---------------------------------------------------------------------------
// planned_action — what a planner decision (StatePlan) forces at one rung. Pure.
// ---------------------------------------------------------------------------

/// The action `plan` gives `node` at `thresh`. A rung this run has already
/// scheduled is Skip: the session's node owns completion. Otherwise the
/// action names the source the metadata holds (best_usable_fd_source_key,
/// best_usable_es_source_key), for the scheduler trace.
inline NodeAction
planned_action(const CalcNode &node, const nlohmann::json &meta, double thresh,
               const StatePlan &plan, const SolvedRungs &solved) {
  if (plan.kind == StatePlan::Kind::load) return NodeAction::Skip;
  const double start = plan.start_thresh;
  const bool   at_start =
      std::abs(thresh - start) <= 1e-12 * std::max(thresh, start);
  if (!at_start && thresh > start) return NodeAction::Skip;
  const int         k   = default_k_for_thresh(thresh);
  const std::string key = protocol_key(thresh, k);
  if (solved.contains(node.id + "|" + key)) return NodeAction::Skip;
  if (at_start && plan.kind == StatePlan::Kind::fresh) return NodeAction::Fresh;
  // An excited-state bundle (D-025): the source the bundle loader picks.
  if (node.kind == CalcKind::ES) {
    const std::string source = best_usable_es_source_key(
        meta, node.tda ? "tda" : "full", "closed_shell", node.n_roots, thresh,
        k, key);
    if (source.empty()) return NodeAction::Fresh;
    return source == key ? NodeAction::Resume : NodeAction::Restart;
  }
  const std::string source = best_usable_fd_source_key(
      meta, node.pert.description(), ResponseMetadata::freq_key(node.freq),
      thresh, k, key);
  if (source.empty()) return NodeAction::Fresh;
  return source == key ? NodeAction::Resume : NodeAction::Restart;
}

// ---------------------------------------------------------------------------
// prerequisites_converged — readiness gate (implicit-by-identity). Pure.
// ---------------------------------------------------------------------------

/// True iff every hard_dep of `node` is converged at protocol step `thresh`. Dep ids
/// are resolved against `dag`; each dep's convergence is checked at the SAME
/// protocol step key (the doc-15 "reconstruct prerequisites' identity at the same
/// protocol_key" contract). A dep id that resolves to no node, or to a node
/// with no converged entry, gates the node out. A dep the planner loaded
/// (`plans`, named path) counts as converged: the planner decided (D-016).
bool prerequisites_converged(const CalcNode &node, const std::vector<CalcNode> &dag,
                const nlohmann::json &meta, double thresh,
                const StatePlans &plans = {}) {
  if (node.prerequisites.empty()) return true;
  const int         k   = default_k_for_thresh(thresh);
  const std::string key = protocol_key(thresh, k);

  auto find_node = [&](const std::string &id) -> const CalcNode * {
    for (const auto &n : dag)
      if (n.id == id) return &n;
    return nullptr;
  };

  for (const auto &dep_id : node.prerequisites) {
    const CalcNode *dep = find_node(dep_id);
    if (!dep) return false;
    const auto plan = plans.find(dep_id);
    if (plan != plans.end() && plan->second.kind == StatePlan::Kind::load &&
        plan->second.met)
      continue;
    bool ok = false;
    if (dep->kind == CalcKind::ES) {
      const auto *e = detail_calc::es_entry(meta, key, *dep);
      ok = e && detail_calc::entry_converged(*e);
    } else {
      const auto *e = detail_calc::fd_entry(meta, dep->pert.description(),
                                            key,
                                            ResponseMetadata::freq_key(dep->freq),
                                            dep->freq);
      ok = e && detail_calc::entry_converged(*e);
    }
    if (!ok) return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// schedule — two-phase, perturbation-parallel waves of WorkItems. Pure.
// ---------------------------------------------------------------------------

/// Emit work in waves (a wave = items with no edges between them, safe to run
/// concurrently). Channels (perturbation / ES bundle) are the unit of
/// parallelism and share every wave. Ordering:
///   Phase A (protocol step = ramp.front()):
///     wave A1 — each perturbation's SEED STATE (ω=0 if present, else lowest |ω|),
///     wave A2 — all remaining nodes at this protocol step.
///   Phase B (protocol step in ramp[1..]): one wave per protocol step of every participating
///     node (embarrassingly parallel).
/// A node participates in a protocol step only if the protocol step is in its protocols ladder
/// and planned_action != Skip. Symbolic nodes (DerivedFD "*") are excluded —
/// they are expansion templates resolved post-ES. Gated nodes (deps not
/// ready) are dropped from the wave and picked up on a later re-schedule
/// (after run() updates the metadata). Empty waves are omitted.
///
/// CONTRACT — only `waves.front()` has FINAL actions. Every wave is reconciled
/// against the SAME metadata snapshot, so a later wave's NodeActions are
/// provisional: e.g. a Phase-B item may read `Fresh` here because P0 has not
/// run yet, when it should be `Restart` once P0 lands. Membership is correct;
/// actions downstream of the front wave are not. This is safe because
/// `CalcManager::run` consumes ONLY the front wave and then re-schedules
/// against fresh metadata — so each wave's actions are final by the time it
/// becomes the front. Do not consume non-front waves' actions directly. (The
/// full list is returned anyway: it is cheap and makes the schedule legible /
/// unit-testable.)
///
/// `plans` (D-016): the session's plan for every schedulable node, which
/// takes planned_action's verdict; a node without a plan is an error.
/// `solved` lists the rungs this run has scheduled.
std::vector<std::vector<WorkItem>>
schedule(const std::vector<CalcNode> &dag,
         const std::vector<double> &global_ramp,
         const nlohmann::json &meta, const StatePlans &plans,
         const SolvedRungs &solved = {}) {
  using detail_calc::ramp_contains;
  std::vector<std::vector<WorkItem>> waves;

  const std::vector<double> ramp = detail_calc::sorted_ramp(global_ramp);
  if (ramp.empty()) return waves;

  // Schedulable (non-symbolic) nodes only.
  std::vector<const CalcNode *> nodes;
  for (const auto &n : dag)
    if (!n.is_symbolic()) nodes.push_back(&n);

  // Does node participate at this protocol step (in ladder, not Skip, deps ready)?
  auto make_item = [&](const CalcNode *n, double thresh,
                       WorkItem &out) -> bool {
    if (!ramp_contains(n->protocols, thresh)) return false;
    const auto plan = plans.find(n->id);
    if (plan == plans.end())
      throw std::logic_error("response state " + n->id + " has no plan");
    if (!prerequisites_converged(*n, dag, meta, thresh, plans)) return false;
    const NodeAction a = planned_action(*n, meta, thresh, plan->second, solved);
    if (a == NodeAction::Skip) return false;
    out = WorkItem{n, thresh, a};
    return true;
  };

  // ---- Phase A: seed states then the rest, at ramp.front() --------------------
  const double p0 = ramp.front();

  // Per-perturbation seed state = the node with freq closest to 0 (static
  // preferred), among nodes that participate at p0. Precomputed once (O(n));
  // std::map keeps A1 in a deterministic (rank-consistent) perturbation order.
  std::map<std::string, const CalcNode *> seed_by_perturbation;
  for (const CalcNode *n : nodes) {
    if (!ramp_contains(n->protocols, p0)) continue;
    const std::string ch = perturbation_key(*n);
    auto it = seed_by_perturbation.find(ch);
    if (it == seed_by_perturbation.end() ||
        std::abs(n->freq) < std::abs(it->second->freq))
      seed_by_perturbation[ch] = n;
  }

  std::vector<WorkItem> waveA1, waveA2;
  for (const auto &[ch, seed] : seed_by_perturbation) {
    WorkItem it;
    if (make_item(seed, p0, it)) waveA1.push_back(it);
  }
  // A2: every other p0 node that isn't its perturbation's seed state.
  for (const CalcNode *n : nodes) {
    auto a = seed_by_perturbation.find(perturbation_key(*n));
    if (a != seed_by_perturbation.end() && a->second == n) continue;
    WorkItem it;
    if (make_item(n, p0, it)) waveA2.push_back(it);
  }
  if (!waveA1.empty()) waves.push_back(std::move(waveA1));
  if (!waveA2.empty()) waves.push_back(std::move(waveA2));

  // ---- Phase B: one wave per finer protocol step -----------------------------------
  for (size_t i = 1; i < ramp.size(); ++i) {
    std::vector<WorkItem> wave;
    for (const CalcNode *n : nodes) {
      WorkItem it;
      if (make_item(n, ramp[i], it)) wave.push_back(it);
    }
    if (!wave.empty()) waves.push_back(std::move(wave));
  }

  return waves;
}

// ---------------------------------------------------------------------------
// fd_leg_accuracy — the honest verdict of one FD leg behind a property row. Pure.
// ---------------------------------------------------------------------------

/// Accuracy record of the FD state (pert, freq) a property row was built from,
/// at target (thresh, k): the SAME source the loader uses (best_usable_fd_source_key),
/// so the reported accuracy is the accuracy of the state actually contracted.
/// `converged` is the STRICT verdict — an accepted-at-maxiter state (converged
/// forced true in fd_states so dependent gates unblock) reports converged=false
/// and accepted=true. Shared by the alpha, beta/Raman and 2PA assemblers
/// (review findings C9, D4). Returns {input, freq, source_protocol_key,
/// converged, accepted, bsh_residual}; source_protocol_key is "" when no usable
/// state exists. A leg the planner loaded (`plans`, named path, D-016) adds
/// met = true and decided_by = "planner"; `converged` stays the engine's
/// recorded verdict, so a reader can tell the two apart (B8).
inline nlohmann::json fd_leg_accuracy(const nlohmann::json &meta, const Perturbation &p,
                                      double freq, double thresh, int k,
                                      const std::string &key,
                                      const StatePlans &plans = {}) {
  const std::string chan = p.description();
  const std::string fk   = ResponseMetadata::freq_key(freq);
  const std::string sk   = best_usable_fd_source_key(meta, chan, fk, thresh, k, key);
  nlohmann::json a{{"input", chan}, {"freq", freq}, {"source_protocol_key", sk},
                   {"converged", false}, {"accepted", false}, {"bsh_residual", 0.0}};
  if (!sk.empty()) {
    const auto &e = meta["fd_states"][chan][sk][fk];
    const bool acc    = e.value("accepted", false);
    a["converged"]    = e.value("converged", false) && !acc;
    a["accepted"]     = acc;
    a["bsh_residual"] = e.value("bsh_residual", 0.0);
    // Only the state the planner judged: a loaded tighter rung read through
    // a coarser source here is not it.
    const auto plan = plans.find(fd_node_id(p, freq));
    if (plan != plans.end() && plan->second.kind == StatePlan::Kind::load &&
        plan->second.met && plan->second.loaded_key == sk) {
      a["met"]        = true;
      a["decided_by"] = "planner";
    }
  }
  return a;
}

// ---------------------------------------------------------------------------
// audit_dropped_work — planned work the run ended without delivering. Pure.
// ---------------------------------------------------------------------------

/// An empty schedule does NOT mean every planned node executed: honest-climb can
/// walk a stubborn prerequisite up the whole ladder unconverged, and nodes gated
/// on it are then invisible to schedule(). Audit the plan against the metadata.
/// Returns [{id, top_protocol_key, reason}]. The run records this under
/// run_summary/dropped_work and names it in stop_reason.
inline nlohmann::json audit_dropped_work(const std::vector<CalcNode> &dag,
                                         const nlohmann::json &meta,
                                         const std::set<std::string> &stalled_ids) {
  nlohmann::json dropped = nlohmann::json::array();
  for (const auto &n : dag) {
    if (n.kind != CalcKind::VBC || n.protocols.empty()) continue;
    const std::string top = protocol_key_at(n.protocols.back());
    const bool has_entry = meta.contains("vbc_states") &&
                           meta["vbc_states"].contains(n.id) &&
                           meta["vbc_states"][n.id].contains(top);
    if (has_entry && meta["vbc_states"][n.id][top].value("converged", false))
      continue;
    // solve_vbc is one-shot and always saves converged=true, so a present VBC
    // entry is caught by the `continue` above; the "built ... not converged" arm
    // is defensive coverage in case VBC ever gains partial saves.
    const char *reason =
        stalled_ids.count(n.id)
            ? "stalled (quarantined by the no-progress guard)"
            : has_entry ? "built at the top rung but not converged"
                        : "prerequisites never converged (gated out of every wave)";
    dropped.push_back({{"id", n.id}, {"top_protocol_key", top}, {"reason", reason}});
  }
  // Two-photon legs ("*" DerivedFD) expand only once their ES bundle converges
  // at its top rung. If it never did (stalled, budget exhausted, quarantined),
  // the legs and the 2PA rows they feed were never produced: report it.
  for (const auto &n : dag) {
    if (!n.is_symbolic() || n.prerequisites.empty()) continue;
    const CalcNode *es = nullptr;
    for (const auto &e : dag)
      if (e.id == n.prerequisites.front()) { es = &e; break; }
    if (!es || es->protocols.empty()) continue;
    const std::string top = protocol_key_at(es->protocols.back());
    const auto *entry = detail_calc::es_entry(meta, top, *es);
    if (entry && detail_calc::entry_converged(*entry)) continue;
    dropped.push_back({{"id", n.id}, {"top_protocol_key", top},
                       {"reason", "excited-state bundle " + es->id +
                                      " never converged at its top rung "
                                      "(two-photon legs never expanded)"}});
  }
  return dropped;
}

// ---------------------------------------------------------------------------
// Executor seam (World-bound impl + CalcManager::run land in the next
// increment, calc/calc_executor.hpp). Declared here so the scheduling side
// and its tests agree on the contract.
// ---------------------------------------------------------------------------

struct NodeResult {
  bool                converged = false;
  std::string         reached_protocol_key;
  std::vector<double> es_root_freqs;   // vestigial: DerivedFD expansion is now
                                       // metadata-driven (expand_converged_es
                                       // reads roots from response_metadata.json).
                                       // run()'s loop discards NodeResult, so this
                                       // is no longer populated or read.
};

/// Runs one protocol step of one node: load nearest seed -> iterate at this protocol
/// -> save (+metrics). The only World-bound surface of the calc manager.
struct ICalcExecutor {
  virtual NodeResult run_protocol(const WorkItem &item) = 0;
  virtual ~ICalcExecutor() = default;
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_CALC_CALC_MANAGER_HPP
