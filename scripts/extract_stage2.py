#!/usr/bin/env python3
"""
Reads the six Stage 2 run folders and fills results/stage2_results.csv.

    python scripts/extract_stage2.py <folder holding the runs>

The folder is the one with stage2_a30_mu000, stage2_a30_mu010 ... inside it.
Each of those needs spcforc, nodout and glstat; anything missing is reported and
the row is left blank rather than guessed.

Nothing here is a new measurement. It reads the same three files LS-PrePost
reads and takes the same numbers off them, which is the point: the extraction
is written down instead of clicked, so the numbers in the CSV can be traced
back to a line in a file.

One thing worth knowing about spcforc: every output block already ends with

         force resultants   =   -2.9730E-02   -3.7637E-14    3.7950E+00

which is the x, y, z sum over all 35 constrained root nodes. So W and P are in
the file directly, one line per output time, and there is no summing to do.
"""

import csv
import math
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

# The six runs, in the order they go into the CSV.
CASES = [
    ("a30_mu000", 30.0, 0.00),
    ("a30_mu010", 30.0, 0.10),
    ("a30_mu020", 30.0, 0.20),
    ("a30_mu030", 30.0, 0.30),
    ("a45_mu000", 45.0, 0.00),
    ("a45_mu020", 45.0, 0.20),
]

T_RAMP = 5.0e-3          # end of the ramp
T_END = 5.5e-3           # end of the analysis
TOL = 1.0e-6             # output times land on the nearest cycle, not on the grid

# KE/IE is only read over the second half of the ramp. Earlier than that the
# lance has almost no internal energy - at 1 ms this model holds 2.7e-3 Nmm
# against 1.15 at the end - and the ratio is then set by its own denominator,
# not by anything dynamic. It reads 7.3 % at 1.02 ms and 0.32 % at 2.51 ms on
# the same run, with the kinetic energy itself still rising between the two.
# The absolute peak and its time are printed for every run so the early part is
# still visible.
T_KE_FROM = 2.5e-3

FLOAT = re.compile(r"[-+]?\d*\.\d+[Ee][-+]\d+")


def read_spcforc(path):
    """
    Returns [(time, fx, fy, fz)] - one entry per output block, from the
    'force resultants' line that closes it.
    """
    out = []
    t = None
    with open(path, errors="ignore") as f:
        for line in f:
            if "output at time" in line:
                m = FLOAT.search(line)
                t = float(m.group()) if m else None
            elif "force resultants" in line and t is not None:
                v = FLOAT.findall(line)
                if len(v) >= 3:
                    out.append((t, float(v[0]), float(v[1]), float(v[2])))
    return out


def read_nodout(path):
    """
    Returns (tip_node_id, [(time, z_disp)]) for the node furthest along x.

    The tip node is found rather than hard-coded: in the first block the record
    with the largest x coordinate is the free end. That survives a change of
    mesh density, which a hard-coded id would not.
    """
    times, tip, first = [], None, {}
    t, t0 = None, None
    in_block = False
    with open(path, errors="ignore") as f:
        for line in f:
            if "n o d a l" in line:
                if tip is None and first:
                    tip = max(first, key=lambda n: first[n][0])
                    times.append((t0, first[tip][1]))
                m = re.search(r"at time\s+([-+]?\d*\.\d+[Ee][-+]\d+)", line)
                t = float(m.group(1)) if m else None
                if t0 is None:
                    t0 = t
                in_block = False
            elif "nodal point" in line:
                in_block = True
            elif in_block:
                m = re.match(r"\s*(\d+)\s", line)
                if not m:
                    in_block = False
                    continue
                nid = int(m.group(1))
                v = [float(x) for x in FLOAT.findall(line)]
                if len(v) < 12:
                    continue
                if tip is None:
                    first[nid] = (v[9], v[2])     # x coordinate, z displacement
                elif nid == tip:
                    times.append((t, v[2]))
    if tip is None and first:                     # a file with a single block
        tip = max(first, key=lambda n: first[n][0])
        times.append((t0, first[tip][1]))
    return tip, times


def read_glstat(path):
    """Returns [{label: value}] - one dict per output block."""
    blocks, cur = [], {}
    label = re.compile(r"\s*([a-z].*?)\.{2,}\s*([-+]?\d*\.\d+[Ee][-+]\d+)\s*$")
    with open(path, errors="ignore") as f:
        for line in f:
            if line.strip().startswith("dt of cycle"):
                if cur:
                    blocks.append(cur)
                cur = {}
                continue
            m = label.match(line)
            if m:
                cur[m.group(1).strip()] = float(m.group(2))
    if cur:
        blocks.append(cur)
    return [b for b in blocks if "time" in b]


def wedge_ratio(angle_deg, mu):
    """Closed form, repeated here so the windows can be printed with a target."""
    ta = math.tan(math.radians(angle_deg))
    denom = 1.0 - mu * ta
    return float("inf") if denom <= 1e-12 else (ta + mu) / denom


def ratio_over(forces, t0, t1):
    """
    Mean W, mean P and W/P over a time window. Ratio of the means, not the mean
    of the ratios - see postprocess_stage2.measured_ratio.
    """
    w = [abs(r[1]) for r in forces if t0 - TOL <= r[0] <= t1 + TOL]
    p = [abs(r[3]) for r in forces if t0 - TOL <= r[0] <= t1 + TOL]
    if not w or mean(p) == 0.0:
        return 0, float("nan"), float("nan"), float("nan")
    return len(w), mean(w), mean(p), mean(w) / mean(p)


# The window W/P is measured in, and the reason it is not the hold.
#
# W/P = (tan a + mu)/(1 - mu tan a) is a sliding relation: it comes from Coulomb
# friction with relative motion at the interface. During the hold the plate is
# stationary, so there is no relative motion and the contact is stuck. A stuck
# penalty interface carries whatever tangential force its stick spring was
# holding when motion stopped, and sheds some of it as the lance settles back.
# Measured in the hold the four friction runs come out 0.6 to 1.9 % low, with
# the miss growing with mu and both frictionless runs unaffected - which is the
# signature of a tangential force relaxing, since with mu = 0 there is none to
# relax. Measured while the plate is moving all six sit on the closed form.
#
# The plate follows a quintic smoothstep, so its speed peaks at 2.5 ms and is
# below a tenth of that by 4.5 ms. 2.5 to 3.5 ms is the fastest part of the
# stroke, with the force already at half its final value. The choice is not
# delicate: over all 35 windows between 0.8 and 1.5 ms long inside 2.2 to
# 4.2 ms, every case has a median deviation inside 0.2 % and no window in any
# case exceeds 0.93 %.
T_SLIDE = (2.5e-3, 3.5e-3)

# Also printed for every run, so the sliding and stationary readings stay side
# by side and the difference above can be seen rather than taken on trust.
WINDOWS = [
    ("2.0-3.0 ms", 2.0e-3, 3.0e-3),
    ("2.5-3.5 ms", T_SLIDE[0], T_SLIDE[1]),
    ("3.0-4.0 ms", 3.0e-3, 4.0e-3),
    ("4.0-5.0 ms", 4.0e-3, 5.0e-3),
    ("hold",       T_RAMP, T_END),
]


def nearest(series, t, key=lambda r: r[0], window=5.0e-5):
    """
    Closest entry to time t, or None if nothing lands near it. The window is
    five output intervals: an incomplete run then leaves the cell blank instead
    of quietly reporting whatever it managed to write before it stopped.
    """
    if not series:
        return None
    best = min(series, key=lambda r: abs(key(r) - t))
    return best if abs(key(best) - t) <= window else None


def mean(v):
    return sum(v) / len(v)


def stdev(v):
    if len(v) < 2:
        return float("nan")
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def find_run(root, key):
    """The run folder for one case, by name, one level down or two."""
    for d in sorted(root.rglob("*")):
        if d.is_dir() and key in d.name and (d / "spcforc").exists():
            return d
    return None


def extract(folder):
    """The three files of one run, reduced to the row that goes in the CSV."""
    row, notes = {}, []

    forces = read_spcforc(folder / "spcforc")
    row["_history"] = forces
    row["_windows"] = [(name, *ratio_over(forces, a, b)) for name, a, b in WINDOWS]
    slide = [r for r in forces if T_SLIDE[0] - TOL <= r[0] <= T_SLIDE[1] + TOL]
    if not slide:
        notes.append("no spcforc output inside the sliding window")
    else:
        w = [abs(r[1]) for r in slide]
        p = [abs(r[3]) for r in slide]
        row["W_N"] = mean(w)
        row["P_N"] = mean(p)
        # Reported, not used: the comparison divides the means, which is the
        # impulse ratio over the window. Dividing sample by sample and averaging
        # that biases the answer - see postprocess_stage2.measured_ratio.
        per_sample = [a / b for a, b in zip(w, p) if b > 0.0]
        if per_sample:
            row["ratio_scatter_pct"] = stdev(per_sample) / mean(per_sample) * 100.0
        notes.append(f"{len(slide)} samples while sliding, "
                     f"{slide[0][0] * 1e3:.3f} to {slide[-1][0] * 1e3:.3f} ms")

    # The hold is kept because it is the end of the stroke: the forces there are
    # the fully developed ones, which is what the frictional work is normalised
    # on. It is not what W/P is compared on - see the note above T_SLIDE.
    hold = [r for r in forces if T_RAMP - TOL <= r[0] <= T_END + TOL]
    if hold:
        row["W_hold_N"] = mean([abs(r[1]) for r in hold])
        row["P_hold_N"] = mean([abs(r[3]) for r in hold])
    else:
        notes.append("no spcforc output inside the hold window")

    if (folder / "nodout").exists():
        tip, disp = read_nodout(folder / "nodout")
        at_ramp = nearest(disp, T_RAMP)
        if at_ramp:
            row["tip_disp_mm"] = at_ramp[1]
            notes.append(f"tip node {tip}, z at {at_ramp[0] * 1e3:.3f} ms")
    else:
        notes.append("no nodout")

    if (folder / "glstat").exists():
        g = read_glstat(folder / "glstat")
        window = [b for b in g if T_KE_FROM - TOL <= b["time"] <= T_RAMP + TOL]
        ratios = [b["kinetic energy"] / b["internal energy"] * 100.0
                  for b in window if b.get("internal energy", 0.0) > 0.0]
        if ratios:
            row["ke_over_ie_pct"] = max(ratios)
        ramp = [b for b in g if b["time"] <= T_RAMP + TOL]
        if ramp:
            pk = max(ramp, key=lambda b: b.get("kinetic energy", 0.0))
            notes.append(f"peak kinetic energy {pk['kinetic energy']:.3e} Nmm "
                         f"at {pk['time'] * 1e3:.3f} ms")
        end = nearest(g, T_RAMP, key=lambda b: b["time"])
        if end and end.get("internal energy", 0.0) > 0.0:
            row["internal_energy_Nmm"] = end["internal energy"]
            row["sliding_over_ie_pct"] = (end.get("sliding interface energy", 0.0)
                                          / end["internal energy"] * 100.0)
        if end:
            row["energy_ratio"] = end.get("total energy / initial energy", float("nan"))
        if g:
            worst = max(g, key=lambda b: abs(b.get("total energy / initial energy", 1.0) - 1.0))
            notes.append("energy ratio worst "
                         f"{worst.get('total energy / initial energy', float('nan')):.5f} "
                         f"at {worst['time'] * 1e3:.3f} ms")
        else:
            notes.append("glstat had no output blocks")
    else:
        notes.append("no glstat")

    return row, notes


HEADER = """\
# Stage 2 results, one row per run. Written by scripts/extract_stage2.py - do
# not edit by hand, rerun the script instead.
#
# tip_disp_mm           tip node z displacement at the end of the ramp, from
#                       nodout. It is close to 0.600 but not equal to it and can
#                       land either side: the plate face is inclined, so the
#                       face height under the tip depends on where the tip
#                       corner has moved to, and the corner both draws in and
#                       rides forward as the lance bends. Recorded, not assumed.
# W_N                   summed x reaction at the root, spcforc 'force
#                       resultants', mean over 2.50 to 3.50 ms - the fastest
#                       part of the stroke, where the contact is sliding. W/P is
#                       a sliding relation and this is where it is read.
# P_N                   summed z reaction at the root, same line, same window.
# ratio_scatter_pct     standard deviation of the per-sample W/P over the same
#                       window, in percent. For the report only - the comparison
#                       uses W_N/P_N, which is the unbiased estimator.
# W_hold_N, P_hold_N    the same two forces over the hold, 5.00 to 5.50 ms,
#                       where the plate has stopped and the contact is stuck.
#                       These are the fully developed end-of-stroke forces and
#                       the frictional work is normalised on them. W/P read here
#                       comes out 0.6 to 1.9 % low for mu > 0 and unaffected for
#                       mu = 0, which is the stick relaxing, not the friction
#                       law - see the note above T_SLIDE in the script.
# ke_over_ie_pct        glstat, max kinetic/internal energy over 2.5 to 5.0 ms,
#                       the half of the ramp where the internal energy is
#                       developed enough for the ratio to mean something.
# sliding_over_ie_pct   glstat, sliding interface energy / internal energy at
#                       the end of the ramp. With friction on most of this is
#                       real frictional work, so it is not a defect screen -
#                       postprocess_stage2.py compares it against the closed
#                       form instead.
# internal_energy_Nmm   glstat, internal energy at the end of the ramp. Carried
#                       so the sliding energy can be compared in N mm.
# energy_ratio          glstat, total energy / initial energy at the end of the
#                       ramp. This is the check that holds for every case.
#
# Compared against W/P = (tan a + mu)/(1 - mu tan a). The ratio holds no
# stiffness, so the Stage 1 deviation from beam theory does not enter it.
"""

COLUMNS = ["case", "angle_deg", "mu", "tip_disp_mm", "W_N", "P_N",
           "ratio_scatter_pct", "W_hold_N", "P_hold_N", "ke_over_ie_pct",
           "sliding_over_ie_pct", "internal_energy_Nmm", "energy_ratio"]


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "ls-dyna" / "stage2"
    if not root.is_dir():
        raise SystemExit(f"{root} is not a folder.")

    print(f"Reading runs from {root}\n")
    rows, history = [], []
    for key, angle, mu in CASES:
        folder = find_run(root, key)
        if folder is None:
            print(f"{key:>12}  not found")
            rows.append({"case": key, "angle_deg": f"{angle:.0f}", "mu": f"{mu:.2f}"})
            continue

        row, notes = extract(folder)
        print(f"{key:>12}  {folder.name}")
        for n in notes:
            print(f"               {n}")
        if row.get("P_N"):
            print(f"               W {row['W_N']:.4f} N   P {row['P_N']:.4f} N   "
                  f"W/P {row['W_N'] / row['P_N']:.4f}")

        target = wedge_ratio(angle, mu)
        if row.get("_windows"):
            print(f"               {'window':<12}{'n':>4}{'W':>9}{'P':>9}"
                  f"{'W/P':>9}{'vs ' + format(target, '.4f'):>12}")
            for name, n, w, p, ratio in row["_windows"]:
                if not n:
                    continue
                dev = (ratio / target - 1.0) * 100.0 if target else float("nan")
                print(f"               {name:<12}{n:>4}{w:>9.4f}{p:>9.4f}"
                      f"{ratio:>9.4f}{dev:>11.2f}%")
        if row.get("_history"):
            history.append((key, row.pop("_history")))
        row.pop("_windows", None)
        row.pop("_history", None)
        print()

        out = {"case": key, "angle_deg": f"{angle:.0f}", "mu": f"{mu:.2f}"}
        fmt = {"tip_disp_mm": "{:.4f}", "W_N": "{:.4f}", "P_N": "{:.4f}",
               "ratio_scatter_pct": "{:.1f}", "W_hold_N": "{:.4f}",
               "P_hold_N": "{:.4f}", "ke_over_ie_pct": "{:.4f}",
               "sliding_over_ie_pct": "{:.2f}", "internal_energy_Nmm": "{:.4f}",
               "energy_ratio": "{:.5f}"}
        for c, f in fmt.items():
            out[c] = f.format(row[c]) if c in row else ""
        rows.append(out)

    path = RESULTS / "stage2_results.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        f.write(HEADER)
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"Written  {path}")

    # The whole reaction history, every output time, for all six runs. It costs
    # nothing to write and it is the only way to see whether the ratio is flat
    # while the plate slides or only settles once it stops.
    if history:
        hpath = RESULTS / "stage2_ratio_history.csv"
        with open(hpath, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["case", "time_ms", "W_N", "P_N"])
            for case, series in history:
                for t, fx, _fy, fz in series:
                    w.writerow([case, f"{t * 1e3:.5f}",
                                f"{abs(fx):.6f}", f"{abs(fz):.6f}"])
        print(f"Written  {hpath}")

    done = sum(1 for r in rows if r.get("W_N"))
    print(f"{done} of {len(CASES)} runs extracted. "
          "Next: python scripts/postprocess_stage2.py")


if __name__ == "__main__":
    main()
