#!/usr/bin/env python3
"""
IHIS Phase 7, T4.4 / T8.1-T8.3: blunted cone post-processing.

Usage: bluntedCone_post.py <run dir> [time]

Surface (Casseau et al. 2016, Fig. 2d-f), x = axial distance from the
stagnation point [cm]:
    Cp = (p - p_inf)/q_inf,  Cf = |tau_t|/q_inf,  St = q/(0.5 rho_inf U_inf^3)
    q_inf = 0.5 rho_inf U_inf^2.  St from q_tot and from q_cond (Q17).
References: hy2Foam_*.csv, monacoDSMC_*.csv in references/ next to this script
(filled once by './Allrun import').

Stagnation line (Fig. 2a-c): T/T_inf, rho/rho_inf, U/U_inf against the
distance from the wall [mm]; nexaFoam overlay if
references/nexaFoam/postProcessing/stagLine/0.0002/axis*.csv exists.
T_v is drawn only where it varies (not in 1T, not for a frozen T_v).
The outlet-corner faces (x > XMAX_CM) are left out of the surface plots.

Writes <run dir>/post_surface.png, post_stagLine.png, post_surface.csv.
"""
import csv, glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, UINF, PINF, TINF = 5.113e-4, 2764.5, 21.9139, 144.4
QINF = 0.5*RHO*UINF**2
XMAX_CM = 4.15     # last cone faces at the outlet corner are a post-processing artefact

run = sys.argv[1]
refDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "references")
times = sorted((d for d in os.listdir(run) if re.fullmatch(r"[0-9.eE+-]+", d) and d != "0"),
               key=float)
t = sys.argv[2] if len(sys.argv) > 2 else times[-1]
print(f"run {run}, time {t}")

def patchBlock(fname, patch="cone"):
    s = open(f"{run}/{t}/{fname}").read()
    i = s.index("boundaryField")
    j = re.search(r"\n\s*" + patch + r"\s*\n\s*\{", s[i:]).end() + i
    depth, k = 1, j
    while depth:
        depth += {"{": 1, "}": -1}.get(s[k], 0)
        k += 1
    return s[j:k]

def valueList(block, vector=False):
    m = re.search(r"\bvalue\s+nonuniform\s+List<\w+>\s*\d+\s*\((.*?)\)\s*;", block, re.S)
    body = m.group(1)
    if vector:
        return np.array([[float(v) for v in g.split()] for g in re.findall(r"\(([^()]*)\)", body)])
    return np.array([float(v) for v in body.split()])

# ---- surface ---------------------------------------------------------------
Cf_ = valueList(patchBlock("C"), vector=True)
p = valueList(patchBlock("p"))
xs = Cf_[:, 0].min()

def readCsv(path):
    rows = [r for r in csv.reader(l for l in open(path) if not l.startswith("#"))]
    head, data = rows[0], np.array(rows[1:], dtype=float)
    return {h: data[:, i] for i, h in enumerate(head)}

hf = readCsv(f"{run}/postProcessing/wallHeatFlux/{t}/cone.csv")
ws = readCsv(f"{run}/postProcessing/wallShearStress/{t}/cone.csv")

o = np.argsort(Cf_[:, 0]); xp, cp = (Cf_[o, 0] - xs)*100, (p[o] - PINF)/QINF
o = np.argsort(hf["x"]); xh = (hf["x"][o] - xs)*100
st_tot = hf["q_tot"][o]/(0.5*RHO*UINF**3)
st_cond = hf["q_cond"][o]/(0.5*RHO*UINF**3)
o = np.argsort(ws["x"]); xw, cf = (ws["x"][o] - xs)*100, ws["tau_t"][o]/QINF

with open(f"{run}/post_surface.csv", "w") as f:
    f.write("# x_cm,Cp | x_cm,St_tot,St_cond | x_cm,Cf  (separate sorted columns)\n")
    for a in zip(xp, cp, xh, st_tot, st_cond, xw, cf):
        f.write(",".join(f"{v:.8g}" for v in a) + "\n")

def ref(name):
    path = f"{refDir}/{name}.csv"
    return np.loadtxt(path, delimiter=",") if os.path.exists(path) else None

fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for a, (key, x, y, lab) in zip(ax, [("cp", xp, cp, "Pressure coefficient"),
                                    ("cf", xw, cf, "Friction coefficient"),
                                    ("St", xh, st_tot, "Stanton number")]):
    for src, mk in (("hy2Foam", "s"), ("monacoDSMC", "o")):
        r = ref(f"{src}_{key}")
        if r is not None:
            a.plot(r[:, 0], r[:, 1], mk, mfc="none", ms=5, label=src)
    m = x < XMAX_CM
    a.plot(x[m], y[m], "-", lw=1.6, label=f"IHIS ({os.path.basename(run)})")
    if key == "St":
        m = xh < XMAX_CM
        a.plot(xh[m], st_cond[m], "--", lw=1.0, label="IHIS, conduction only")
    a.set_xlabel("Axial distance from stagnation point [cm]"); a.set_ylabel(lab)
    a.set_xlim(0, 4.5); a.grid(alpha=0.3); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{run}/post_surface.png", dpi=150)
print(f"  wrote {run}/post_surface.png, post_surface.csv")

# ---- stagnation line ---------------------------------------------------------
def stagLine(case, tt):
    files = glob.glob(f"{case}/postProcessing/stagLine/{tt}/axis*.csv")
    if not files:
        return None
    d = {}
    for f in files:
        d.update(readCsv(f))
    return d

xs_mm = xs*1000
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for case, tt, lab, ls in ((run, t, "IHIS", "-"), (f"{refDir}/nexaFoam", "0.0002", "nexaFoam", "--")):
    d = stagLine(case, tt)
    if d is None:
        print(f"  no stagLine output for {lab} ({case}, {tt})")
        continue
    xx = d["x"]*1000 - xs_mm
    ucols = [k for k in d if re.fullmatch(r"U_?[0-2xyzXYZ]", k)]
    Umag = np.sqrt(sum(d[k]**2 for k in ucols)) if ucols else np.zeros_like(xx)
    if not ucols:
        print(f"  {lab}: no U columns in {sorted(d)}")
    ax[0].plot(xx, d["TTR"]/TINF, ls, label=f"{lab} T_tr")
    if "TVib" in d and np.ptp(d["TVib"]) > 1e-6*TINF:   # skip 1T / frozen T_v
        ax[0].plot(xx, d["TVib"]/TINF, ls, lw=0.8, label=f"{lab} T_v")
    ax[1].plot(xx, d["rho"]/RHO, ls, label=lab)
    ax[2].plot(xx, Umag/UINF, ls, label=lab)
for a, lab in zip(ax, ["T / T_inf", "rho / rho_inf", "U / U_inf"]):
    a.set_xlabel("Stagnation line position [mm]"); a.set_ylabel(lab)
    a.set_xlim(-3, 0); a.grid(alpha=0.3); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{run}/post_stagLine.png", dpi=150)
print(f"  wrote {run}/post_stagLine.png")
