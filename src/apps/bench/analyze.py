#!/usr/bin/env python3
"""
analyze.py -- summarize bench_mra_ops JSONL output.

Reads the records emitted by bench_mra_ops and reports, per operation:
  * time vs n at fixed (k, thresh)   -> fitted scaling exponent  time ~ n^p
  * time vs k at fixed (n, thresh)   -> fitted scaling exponent  time ~ k^q
  * communication (messages, MB) and node imbalance per point

Pure stdlib (no numpy). Usage:
    python3 analyze.py bench.jsonl
    python3 analyze.py bench.jsonl --op=apply
"""
import json
import math
import sys
from collections import defaultdict


def loglog_slope(xs, ys):
    """Least-squares slope of log(y) vs log(x); None if <2 usable points."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(pts) < 2:
        return None
    n = len(pts)
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-300:
        return None
    return (n * sxy - sx * sy) / denom


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    if not args:
        print("usage: analyze.py bench.jsonl [--op=NAME]")
        sys.exit(1)

    rows, meta = [], {}
    with open(args[0]) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("record") == "meta":
                meta = r
            else:
                rows.append(r)

    if meta:
        print(f"# run: P={meta.get('P')} threads/rank={meta.get('threads')} "
              f"operator={meta.get('operator')} reps={meta.get('reps')}\n")

    op_filter = opts.get("op")
    ops = sorted({r["op"] for r in rows})
    if op_filter:
        ops = [o for o in ops if o == op_filter]

    for op in ops:
        sel = [r for r in rows if r["op"] == op]
        print(f"=== {op} ===")

        # time vs n, grouped by (k, thresh)
        by_kt = defaultdict(list)
        for r in sel:
            by_kt[(r["k"], r["thresh"])].append(r)
        print("  time vs n:")
        for (k, th), rs in sorted(by_kt.items()):
            rs.sort(key=lambda r: r["n"])
            ns = [r["n"] for r in rs]
            ts = [r["time_s"]["median"] for r in rs]
            p = loglog_slope(ns, ts)
            series = " ".join(f"n={r['n']}:{r['time_s']['median']:.3e}s" for r in rs)
            pstr = f"  ~n^{p:.2f}" if p is not None else ""
            print(f"    k={k} thresh={th:g}: {series}{pstr}")

        # time vs k, grouped by (n, thresh)
        by_nt = defaultdict(list)
        for r in sel:
            by_nt[(r["n"], r["thresh"])].append(r)
        print("  time vs k:")
        for (n, th), rs in sorted(by_nt.items()):
            rs.sort(key=lambda r: r["k"])
            ks = [r["k"] for r in rs]
            ts = [r["time_s"]["median"] for r in rs]
            q = loglog_slope(ks, ts)
            series = " ".join(f"k={r['k']}:{r['time_s']['median']:.3e}s" for r in rs)
            qstr = f"  ~k^{q:.2f}" if q is not None else ""
            print(f"    n={n} thresh={th:g}: {series}{qstr}")

        # communication + imbalance summary
        print("  comm / imbalance:")
        for r in sorted(sel, key=lambda r: (r["k"], r["thresh"], r["n"])):
            mb = r["rmi"]["nbyte_sent"] / 1e6
            print(f"    n={r['n']} k={r['k']} thresh={r['thresh']:g}: "
                  f"msgs={int(r['rmi']['nmsg_sent'])} sent={mb:.1f}MB "
                  f"leaves={int(r['tree']['leaves'])} imbalance={r['tree']['imbalance']:.2f}")
        print()


if __name__ == "__main__":
    main()
