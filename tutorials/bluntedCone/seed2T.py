#!/usr/bin/env python3
"""
Seed the 2T cone start state from a converged 1T TTR field (T7.2):
TVib and TVib_N2 get the 1T TTR internal field (vibration in equilibrium
with translation as the initial guess), with the 2T boundary conditions of
the templates (TVib cone 'calculated' D12, TVib_N2 cone neSmoluchowskiJumpTv).

Usage: seed2T.py <1T TTR file> <TVib.2T template> <TVib_N2.2T template> <0 dir>
"""
import re, sys

ttr, tTv, tTvN2, out = sys.argv[1:5]

s = open(ttr).read()
m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\(.*?\)\s*;", s, re.S)
if not m:
    m = re.search(r"internalField\s+uniform\s+[^;]+;", s)
if not m:
    sys.exit("internalField not found in " + ttr)
internal = m.group(0)

for tmpl, name in ((tTv, "TVib"), (tTvN2, "TVib_N2")):
    t = open(tmpl).read()
    t, n = re.subn(r"internalField\s+uniform\s+\$TValue\s*;", lambda _: internal, t)
    if n != 1:
        sys.exit("template internalField not found in " + tmpl)
    open(f"{out}/{name}", "w").write(t)
    print(f"  wrote {out}/{name} (internalField from {ttr})")
