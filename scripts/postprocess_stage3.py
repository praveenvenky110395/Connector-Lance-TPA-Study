#!/usr/bin/env python3
"""
Turns results/stage3_results.csv into the Stage 3 comparison and figures.

Stage 3 has no closed form of its own. Its prediction is the product of the two
stages already verified - the lance stiffness from Stage 1, the wedge relation
from Stage 2 - and its job is to say what the idealisations cost once the
feature becomes the part.

Sections, in the order they print:

  0  Whether each run's force output may be quoted at all, against the three
     gates recorded in the results file.
  1  Regression - lance_only against Stage 1's FE result times the tooth.
  2  The effective angle at release from the shape, and the force route against
     the shape route station by station.
  3  Insertion against retention, and the asymmetry where both forces pass.
  4  The TPA - blocked or not, and the control run.
  5  Sensitivity - the fine-mesh and slow-ramp runs against the baseline.
"""

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analytical import (Geometry, Material, cycle_forces, effective_angles,
                        tpa_blocking_clearance)
# Deliberately no import from extract_stage3. The gate thresholds come off the
# results file's own header, so what is reported here is what was applied when
# the file was written - and a stale module in a notebook kernel cannot pair a
# new post-processor with an old extractor, which is how this broke once.

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

# Fallbacks only, for a results file written before the limits line existed.
LIMITS = {"min_spc": 8.0, "min_ripple_pct": 5.0, "gap_limit": 1.5,
          "window_spread_limit": 1.5, "windows_us": "50,100,150,250,400"}
_limits = None


def read_limits():
    """The gate thresholds recorded in the results file, if it names them."""
    global _limits
    if _limits is not None:
        return _limits
    _limits = dict(LIMITS)
    path = RESULTS / "stage3_results.csv"
    if path.exists():
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                if "limits:" in line:
                    for part in line.split("limits:", 1)[1].split():
                        if "=" in part:
                            k, v = part.split("=", 1)
                            if k in _limits:
                                _limits[k] = (v if k == "windows_us"
                                              else float(v))
    return _limits

# Same band the other two stages declared. Whether it is met is reported, not
# arranged - Stage 1 missed it and the miss stayed on record.
BAND_PCT = 2.0
ENERGY_BALANCE_PCT = 1.0


def model_spec(case):
    """One case's entry from the index the generator wrote."""
    path = RESULTS / "stage3_model.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())["cases"].get(case)


def load():
    path = RESULTS / "stage3_results.csv"
    if not path.exists():
        raise SystemExit("Run scripts/extract_stage3.py first.")
    num_cols = ("mu", "tpa_gap_mm", "stroke_mm", "tooth_h_mm", "t_ramp_ms",
                "force_dt_us", "samples_per_cycle", "ringing_kHz",
                "ripple_pct", "cutoff_Hz",
                "W_N", "P_N", "W_raw_peak_N", "ring_over_mean", "W_over_P",
                "angle_from_force_deg", "angle_gap_max_deg",
                "angle_gap_mean_deg", "window_spread_deg",
                "angle_from_shape_deg", "rotation_at_release_deg",
                "rotation_at_read_deg", "read_travel_mm", "crest_lift_mm",
                "tip_lift_mm", "release_travel_mm", "W_plateau_N",
                "tpa_engage_travel_mm", "W_at_tpa_engage_N", "ke_over_ie_pct",
                "hourglass_over_ie_pct", "sliding_over_ie_pct", "energy_ratio")
    rows = []
    with open(path) as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            if not r.get("P_N", "").strip():
                continue

            def num(k):
                v = r.get(k, "").strip()
                return float(v) if v else float("nan")

            rows.append({"case": r["case"], "released": r.get("released", ""),
                         "read_at": r.get("read_at", ""),
                         **{k: num(k) for k in num_cols}})
    if not rows:
        raise SystemExit("No Stage 3 results filled in yet.")
    return rows


def force_is_usable(r):
    """
    Whether this run's force columns may be quoted at all.

    Three gates: the output has to resolve the ring, the answer must not move
    with the averaging window, and the two independent routes to the angle must
    agree. The first was set before any run it was applied to. The other two
    were set when the windowed estimator was introduced, after a first look at
    the 0.5 us data, and have not been changed since. A run that fails any of
    them still contributes its displacement measurements, which do not touch
    the force output in any way.
    """
    lim = read_limits()
    spc = r.get("samples_per_cycle", float("nan"))
    ripple = r.get("ripple_pct", float("nan"))
    gap = r.get("angle_gap_max_deg", float("nan"))
    spread = r.get("window_spread_deg", float("nan"))
    if (spc == spc and spc < lim["min_spc"]
            and ripple == ripple and ripple > lim["min_ripple_pct"]):
        return False, f"aliased, {spc:.1f} samples per cycle"
    if spread == spread and spread > lim["window_spread_limit"]:
        return False, f"window dependent, {spread:.1f} deg"
    if gap == gap and gap > lim["gap_limit"]:
        return False, f"force and shape disagree by {gap:.1f} deg"
    return True, ""


def stations():
    """The per-station comparison the extractor wrote, grouped by case."""
    path = RESULTS / "stage3_stations.csv"
    if not path.exists():
        return {}
    out = {}
    with open(path) as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            out.setdefault(r["case"], []).append(
                {"travel": float(r["travel_mm"]),
                 "lift": float(r.get("lift_mm") or 0.0),
                 "W_over_P": float(r["W_over_P"]),
                 "from_force": float(r["angle_from_force_deg"]),
                 "from_shape": float(r["angle_from_shape_deg"]),
                 "gap": float(r["gap_deg"])})
    return out


def plot_cycle(rows, g, m):
    """Force against terminal travel for every run."""
    path = RESULTS / "stage3_curve_history.csv"
    if not path.exists():
        return None
    series = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            series.setdefault(r["case"], []).append(
                (float(r["travel_mm"]), float(r["W_N"]), float(r.get("W_raw_N") or r["W_N"])))

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.4),
                                  gridspec_kw={"width_ratios": [2, 1]})
    colours = {"extract_mu000": "#5a646d", "extract_mu020": "#9d4c19",
               "extract_mu030": "#a03c3c", "extract_tpa_g085": "#2f6f6a",
               "extract_tpa_g010": "#7a5aa0", "insert_mu020": "#3c6ea0",
               "extract_mu020_fine": "#c08a2a", "extract_mu020_slow": "#6b8e23"}
    for r in rows:
        pts = sorted(series.get(r["case"], []))
        if not pts:
            continue
        target = ax2 if r["case"] == "extract_tpa_g010" else ax
        col = colours.get(r["case"], "#333333")
        # The unfiltered trace behind the filtered one, so the size of the
        # ring is visible rather than quietly removed.
        target.plot([p[0] for p in pts], [p[2] for p in pts], "-",
                    color=col, lw=0.5, alpha=0.22)
        target.plot([p[0] for p in pts], [p[1] for p in pts], "-",
                    color=col, lw=1.4,
                    label=r["case"].replace("extract_", "").replace("_", " "))

    ax.set_xlabel("terminal travel  [mm]")
    ax.set_ylabel("pull-out force  W  [N]")
    ax.set_title("Stage 3, the locking cycle")
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(fontsize=8, loc="best", framealpha=0.92)

    ax2.set_xlabel("terminal travel  [mm]")
    ax2.set_title("blocked by the TPA")
    ax2.grid(alpha=0.25, lw=0.6)
    gaps = [r for r in rows if r["case"] == "extract_tpa_g010"]
    if gaps and not math.isnan(gaps[0]["tpa_engage_travel_mm"]):
        ax2.axvline(gaps[0]["tpa_engage_travel_mm"], color="#333333", lw=1.0, ls=":")
        ax2.text(gaps[0]["tpa_engage_travel_mm"], ax2.get_ylim()[1] * 0.5,
                 "  lance meets\n  the TPA", fontsize=7.5, va="center")
    ax2.legend(fontsize=8, loc="best", framealpha=0.92)
    fig.tight_layout()

    out = RESULTS / "figures" / "stage3_cycle.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


def main():
    g, m = Geometry(), Material()
    rows = load()
    a_pred, b_pred = effective_angles(g, m)
    by = {r["case"]: r for r in rows}

    print("Stage 3 - the locking cycle\n")
    print("Prediction is Stage 1 times Stage 2. No new closed form, and none")
    print("needed: what Stage 3 adds is what the idealisations cost.\n")

    # 0 -----------------------------------------------------------------
    print("0. Is the force output usable at all")
    print(f"   {'case':>18} {'dt us':>7} {'samples/cyc':>12} {'ring kHz':>9} "
          f"{'ripple %':>9} {'gap':>7} {'window':>7}  verdict")
    voided = []
    for r in rows:
        ok, why = force_is_usable(r)
        if not ok:
            voided.append(r["case"])
        spc = r["samples_per_cycle"]
        def d(x):
            return f"{x:6.2f}d" if x == x else f"{chr(45):>7}"
        print(f"   {r['case']:>18} {r['force_dt_us']:>7.2f} "
              f"{spc:>12.1f} {r['ringing_kHz']:>9.1f} {r['ripple_pct']:>9.2f} "
              f"{d(r['angle_gap_max_deg'])} {d(r['window_spread_deg'])}  "
              f"{'ok' if ok else why}")
    lim = read_limits()
    ws = lim["windows_us"].split(",")
    print(f"   Gates, as recorded in the results file: at least "
          f"{lim['min_spc']:.0f} samples per cycle of")
    print(f"   the ring where the ripple is over {lim['min_ripple_pct']:.0f} %; "
          f"under {lim['window_spread_limit']:.1f} deg of movement across")
    print(f"   windows of {ws[0]}-{ws[-1]} us; and the force and shape routes "
          f"within {lim['gap_limit']:.1f} deg")
    print("   at every station.")
    if voided:
        print("   Force columns are NOT quoted for: " + ", ".join(voided))
        print("   Their displacement measurements are unaffected and still stand.")

    # 1 -----------------------------------------------------------------
    print("\n1. Regression - does the toothed lance still behave like Stage 1")
    lo, spec = by.get("lance_only"), model_spec("lance_only")
    if lo is None or spec is None:
        print("   lance_only not extracted. Nothing below is trustworthy without")
        print("   it: the tooth is new geometry, and this is the run that says")
        print("   the beam underneath it is still the beam Stage 1 verified.")
    else:
        dev = (lo["P_N"] / spec["P_predicted"] - 1.0) * 100.0
        print(f"   tip lifted {spec['tip_lift']:.3f} mm, nothing else touching the lance")
        print(f"   root reaction   {lo['P_N']:.4f} N")
        print(f"   predicted       {spec['P_predicted']:.4f} N   deviation {dev:+.2f} %")
        print(f"   baseline        {spec.get('baseline_is', '')}")
        print(f"   Stage 1 FE      {spec.get('P_stage1_fe', float('nan')):.4f} N, "
              f"times {spec['tooth_stiffening_pct']:+.2f} % for the tooth")
        print(f"   beam theory     {spec.get('P_beam_theory', float('nan')):.4f} N - "
              "not the baseline, because Stage 1 already")
        print("                   measured and recorded the gap between the 3D FE and")
        print("                   beam theory. Charging it here would count it twice.")

    # 2 -----------------------------------------------------------------
    print("\n2. The effective angle at release, two independent ways")
    print("   The rotation grows with the deflection, so the effective angle is not")
    print("   one number - it starts at the drawn angle and shallows through the")
    print("   stroke. Release is the instant with a sharp definition, so that is")
    print("   where it is quoted.")
    print(f"\n   {'case':>18} {'drawn':>6} {'rotation':>9} {'from shape':>11} "
          f"{'from force':>11} {'predicted':>10} {'shape-pred':>11}")
    for r in rows:
        if math.isnan(r["angle_from_shape_deg"]):
            continue
        ins = r["case"].startswith("insert")
        pred = a_pred if ins else b_pred
        ok, _ = force_is_usable(r)
        fr = (f"{r['angle_from_force_deg']:>10.2f}d"
              if ok and not math.isnan(r["angle_from_force_deg"]) else f"{'void':>11}")
        print(f"   {r['case']:>18} {30.0 if ins else 45.0:>5.0f}d "
              f"{r['rotation_at_release_deg']:>8.2f}d "
              f"{r['angle_from_shape_deg']:>10.2f}d {fr} {pred:>9.2f}d "
              f"{r['angle_from_shape_deg'] - pred:>+10.2f}d")
    print("   The shape column is two node displacements and no force at all, so it")
    print("   does not depend on the filter. That is why it is the primary route.")

    st = stations()
    if st:
        print("\n   Station by station - the whole sweep, both routes. The force column")
        print("   is the ratio of the MEAN axial to the MEAN transverse root reaction")
        print("   over a window, which is Stage 2's estimator. Reading it at a point")
        print("   does not work here: the ripple is most of the force's own level.")
        print("   Stations are placed by the lance's lift, because that is what says")
        print("   which surface is carrying the contact.")
        cases = [c for c in (r["case"] for r in rows) if c in st]
        head = sorted({x["travel"] for c in cases for x in st[c]})
        print(f"\n   {'case':>18} {'lift':>7} {'travel':>8} {'W/P':>8} "
              f"{'force':>9} {'shape':>8} {'gap':>8}")
        for c in cases:
            for x in st[c]:
                print(f"   {c:>18} {x['lift']:>7.3f} {x['travel']:>8.3f} "
                      f"{x['W_over_P']:>8.4f} {x['from_force']:>8.2f}d "
                      f"{x['from_shape']:>7.2f}d {x['gap']:>+7.2f}d")
            print()

    # 3 -----------------------------------------------------------------
    print("\n3. Insertion against retention")
    ins = by.get("insert_mu020")
    ret = by.get("extract_mu020")
    for r in rows:
        if r["case"] == "lance_only" or math.isnan(r["W_N"]):
            continue
        ok, why = force_is_usable(r)
        if not ok:
            print(f"   {r['case']:>18} {'force void':>12}  {why}")
            continue
        if r["released"] == "no" and not math.isnan(r["tpa_gap_mm"]):
            print(f"   {r['case']:>18} {'blocked':>12}  see section 4")
            continue
        w = r["W_plateau_N"] if not math.isnan(r["W_plateau_N"]) else r["W_N"]
        tag = "insertion" if r["case"].startswith("insert") else "retention"
        print(f"   {r['case']:>18} {tag:>12}  W {w:>7.3f} N   "
              f"raw peak {r['W_raw_peak_N']:>7.3f} N  (ring x{r['ring_over_mean']:.1f})")
    if ins and ret and all(force_is_usable(x)[0] for x in (ins, ret)):
        i = ins["W_plateau_N"] if not math.isnan(ins["W_plateau_N"]) else ins["W_N"]
        c = cycle_forces(g, m, ins["mu"])
        print(f"\n   asymmetry measured      {ret['W_N'] / i:.3f}")
        print(f"   asymmetry drawn         {c['asymmetry_drawn']:.3f}   "
              "(hand calculation on the drawn angles)")
        print(f"   asymmetry with rotation {c['asymmetry_rotated']:.3f}   (predicted)")
        print("   If the measured value sits near the rotated prediction and well")
        print("   below the drawn one, the finding holds: the lance's own rotation")
        print("   eats most of the asymmetry the two angles were chosen for.")
    else:
        print("\n   Asymmetry not quoted - it is a ratio of two forces, and at least")
        print("   one of them did not pass section 0.")

    # 4 -----------------------------------------------------------------
    print("\n4. The TPA")
    limit = tpa_blocking_clearance(g, m, g.L)
    print(f"   Blocks below {limit:.3f} mm of clearance - the lance tip lifts")
    print(f"   1.26 times the crest, so the criterion is not the {g.y:.2f} mm")
    print("   protrusion. Clearance is compared against the lift where the TPA is.")
    print(f"\n   {'case':>18} {'gap':>6} {'released':>9} {'engaged at':>11} "
          f"{'W then':>8}   expected")
    for r in rows:
        if math.isnan(r["tpa_gap_mm"]):
            continue
        eng = ("-" if math.isnan(r["tpa_engage_travel_mm"])
               else f"{r['tpa_engage_travel_mm']:.4f}")
        wt = ("-" if math.isnan(r["W_at_tpa_engage_N"])
              else f"{r['W_at_tpa_engage_N']:.3f}")
        print(f"   {r['case']:>18} {r['tpa_gap_mm']:>6.2f} {r['released']:>9} "
              f"{eng:>11} {wt:>8}   "
              f"{'should block' if r['tpa_gap_mm'] < limit else 'should not block'}")
    if "extract_mu020" in by and "extract_tpa_g085" in by:
        a, b = by["extract_mu020"], by["extract_tpa_g085"]
        d = (b["W_N"] / a["W_N"] - 1.0) * 100.0
        dr = (b["release_travel_mm"] / a["release_travel_mm"] - 1.0) * 100.0
        print(f"\n   Control: TPA fitted at 0.85 mm and never touched.")
        print(f"   Release travel differs by {dr:+.3f} %, force by {d:+.3f} %.")
        print("   A contact that is defined but never reached has to cost exactly")
        print("   nothing, and the release travel says so without using the force.")

    # 5 -----------------------------------------------------------------
    base = by.get("extract_mu020")
    var = [r for r in rows if r["case"].startswith("extract_mu020_")]
    if base and var:
        print("\n5. Sensitivity - is this a result or a discretisation?")
        print(f"   {'run':>20} {'tooth h':>8} {'ramp ms':>8} {'release':>9} "
              f"{'vs base':>9} {'angle':>8} {'vs base':>9}")
        print(f"   {'extract_mu020':>20} {base['tooth_h_mm']:>8.3f} "
              f"{base['t_ramp_ms']:>8.1f} {base['release_travel_mm']:>9.4f} {'-':>9} "
              f"{base['angle_from_shape_deg']:>7.2f}d {'-':>9}")
        for r in var:
            dr = (r["release_travel_mm"] / base["release_travel_mm"] - 1.0) * 100.0
            da = r["angle_from_shape_deg"] - base["angle_from_shape_deg"]
            print(f"   {r['case']:>20} {r['tooth_h_mm']:>8.3f} "
                  f"{r['t_ramp_ms']:>8.1f} {r['release_travel_mm']:>9.4f} {dr:>+8.2f}% "
                  f"{r['angle_from_shape_deg']:>7.2f}d {da:>+8.2f}d")
        print("   Release travel and the rotation are the comparison, not the force:")
        print("   they are displacements, so they are the same quantity in all three")
        print("   runs whatever the contact does at the element scale. _fine halves")
        print("   the element over the tooth, _slow doubles the ramp.")

    print("\nSolution quality")
    print(f"   {'case':>18} {'KE/IE %':>9} {'HG/IE %':>9} {'slide/IE %':>11} "
          f"{'energy ratio':>13}")
    for r in rows:
        flag = "" if abs(r["energy_ratio"] - 1.0) * 100 <= ENERGY_BALANCE_PCT \
            else "   <-- energy balance broken"
        print(f"   {r['case']:>18} {r['ke_over_ie_pct']:>9.3f} "
              f"{r['hourglass_over_ie_pct']:>9.3f} {r['sliding_over_ie_pct']:>11.2f} "
              f"{r['energy_ratio']:>13.5f}{flag}")
    print("   Hourglass must be zero - ELFORM -1 is fully integrated. Sliding")
    print("   energy is mostly real frictional work wherever mu > 0, as Stage 2")
    print("   established, so it is reported and not screened.")

    fig = plot_cycle(rows, g, m)
    if fig:
        print(f"\nFigure\n   {fig.relative_to(REPO)}")


if __name__ == "__main__":
    main()
