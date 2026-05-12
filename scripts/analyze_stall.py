#!/usr/bin/env python3
"""analyze_stall.py - per-(protocol, state) stall detector for molresponse runs.

Reads <rundir>/response_base.json, which contains a list of protocols
each with a list of per-iteration BSH residual vectors (one entry per
state). For each (protocol, state) trajectory, find the first iteration
N where the residual hasn't improved by more than (1 - slowdown) over
the past `lookback` iterations.

Usage:
    analyze_stall.py <rundir> [--lookback 3] [--slowdown 0.9] [--no-plot]

Output:
    - Text table to stdout
    - <rundir>/stall_analysis.png unless --no-plot
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def detect_stall(residuals: list[float], lookback: int, slowdown: float) -> int:
    """First iter (index into `residuals`) where progress has stalled.

    "Stalled" means the residual `lookback` iters back was within
    `slowdown` of the current value -- i.e. less than (1 - slowdown)
    relative improvement over that window.

    Returns -1 if the trajectory never stalls within its length.
    """
    for i in range(lookback, len(residuals)):
        if residuals[i] >= slowdown * residuals[i - lookback]:
            return i
    return -1


def load_response_base(rundir: Path) -> dict:
    jf = rundir / "response_base.json"
    if not jf.exists():
        raise FileNotFoundError(jf)
    text = jf.read_text()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        # Some runs write a stray trailing '}'.
        return json.loads(text.rstrip().rstrip("}"))


def trajectories(d: dict) -> list[tuple[float, list[list[float]]]]:
    """Return [(proto_thresh, per_state_trajectories), ...]."""
    out = []
    for proto in d.get("protocol_data", []):
        thresh = proto.get("proto", float("nan"))
        iter_data = proto.get("iter_data", [])
        if not iter_data:
            out.append((thresh, []))
            continue
        n_state = max(
            len(it.get("x_abs_error", {}).get("vals", [])) for it in iter_data
        )
        per_state = [[] for _ in range(n_state)]
        for it in iter_data:
            vals = it.get("x_abs_error", {}).get("vals", [])
            for s in range(n_state):
                if s < len(vals):
                    per_state[s].append(float(vals[s]))
        out.append((thresh, per_state))
    return out


def print_table(rundir: Path, trajs, lookback: int, slowdown: float) -> None:
    print(f"\n{rundir}: {len(trajs)} protocol(s), "
          f"lookback={lookback}, slowdown={slowdown}")
    print(f"  {'proto':>9s} {'state':>5s} {'iters':>5s} "
          f"{'first_res':>11s} {'final_res':>11s} {'best_res':>11s} "
          f"{'stall_iter':>10s}")
    print("  " + "-" * 70)
    for thresh, per_state in trajs:
        if not per_state:
            print(f"  {thresh:>9.0e}   (no iter data)")
            continue
        for s, traj in enumerate(per_state):
            if not traj:
                continue
            stall = detect_stall(traj, lookback, slowdown)
            stall_str = f"{stall}" if stall >= 0 else "-"
            print(f"  {thresh:>9.0e} {s+1:>5d} {len(traj):>5d} "
                  f"{traj[0]:>11.3e} {traj[-1]:>11.3e} {min(traj):>11.3e} "
                  f"{stall_str:>10s}")


def save_plot(rundir: Path, trajs, lookback: int, slowdown: float) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 6))
    colors = plt.get_cmap("tab10").colors

    offset = 0
    proto_lines = []
    legend_added = set()
    for thresh, per_state in trajs:
        if not per_state:
            continue
        n_state = len(per_state)
        proto_iters = max((len(t) for t in per_state), default=0)
        for s, traj in enumerate(per_state):
            if not traj:
                continue
            x = list(range(offset, offset + len(traj)))
            label = f"s{s+1}" if s not in legend_added else None
            legend_added.add(s)
            ax.semilogy(x, traj, marker="o", markersize=4, lw=1.2,
                        color=colors[s % len(colors)], label=label)
            stall = detect_stall(traj, lookback, slowdown)
            if stall >= 0:
                ax.axvline(offset + stall, ls=":",
                           color=colors[s % len(colors)], alpha=0.4)
        proto_lines.append((offset, thresh))
        offset += proto_iters

    # Vertical lines at protocol boundaries (skip the leading 0)
    for o, t in proto_lines[1:]:
        ax.axvline(o, color="k", lw=1.4, alpha=0.6)
    # Label each protocol band at the top.
    ymax = ax.get_ylim()[1]
    for o, t in proto_lines:
        ax.text(o + 0.2, ymax, f"thresh={t:.0e}",
                ha="left", va="top", fontsize=9, alpha=0.7)

    ax.set_xlabel("global iteration (across protocols)")
    ax.set_ylabel("BSH residual (log)")
    ax.set_title(f"{rundir.name}: stall analysis "
                 f"(lookback={lookback}, slowdown={slowdown})")
    ax.legend(loc="upper right")
    out = rundir / "stall_analysis.png"
    plt.tight_layout()
    plt.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rundir", type=Path)
    p.add_argument("--lookback", type=int, default=3,
                   help="compare residual to its value `lookback` iters ago")
    p.add_argument("--slowdown", type=float, default=0.9,
                   help="iter is 'stalled' if r[i] >= slowdown * r[i-lookback]")
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args()

    try:
        d = load_response_base(args.rundir)
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        sys.exit(1)

    trajs = trajectories(d)
    print_table(args.rundir, trajs, args.lookback, args.slowdown)

    if not args.no_plot:
        try:
            out = save_plot(args.rundir, trajs, args.lookback, args.slowdown)
            print(f"\nplot: {out}")
        except ImportError:
            print("\n(matplotlib not available; skipping plot)")


if __name__ == "__main__":
    main()
