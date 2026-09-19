#!/usr/bin/env python3
"""
Reads the Stage 2 results and turns them into the verification numbers.

You fill in results/stage2_results.csv from the LS-DYNA runs, this does the
rest: forms W/P, compares it against the closed-form wedge relation, and draws
the friction figure.

The comparison is on the ratio and not on either force, and the reason is the
whole point of building the study in stages:

    W/P = (tan(angle) + mu) / (1 - mu tan(angle))

holds no stiffness. Not E, not the second moment of area, not the deflection.
So the Stage 1 difference between the FE model and the beam solution - which is
a stiffness difference - cancels out of it, and Stage 2 tests only the contact
normal and the friction law. Stage 1 can be 3.7 % off beam theory and Stage 2
can still be exact, and both statements are true at once.

Where to find each number in the LS-DYNA output
  tip_disp_mm         tip node z displacement at the end of the ramp, from
                      nodout. Magnitude. Close to 0.600 but not equal to it, and
                      it can land either side: the face is inclined, so the
                      height of the face under the tip depends on where the tip
                      corner has got to, and that corner draws in and rides
                      forward as the lance bends. Recorded, not assumed.
  W_N                 summed x reaction at the root, spcforc, averaged over
                      2.50 to 3.50 ms - the fastest part of the stroke, where
                      the contact is sliding. See measured_ratio().
  P_N                 summed z reaction at the root, same file, same window.
  ratio_scatter_pct   standard deviation of the per-sample W/P over the same
                      window, as a percent. Recorded for the report; the
                      comparison uses W_N/P_N.
  W_hold_N, P_hold_N  the same two forces over the hold, where the plate has
                      stopped. End-of-stroke values; the frictional work is
                      normalised on them and they are reported alongside.
  ke_over_ie_pct      glstat, max kinetic/internal energy over 2.5 to 5.0 ms.
  sliding_over_ie_pct glstat, sliding interface energy over internal energy at
                      the end of the ramp.
  energy_ratio        glstat, total energy over initial energy at the same
                      instant.

scripts/extract_stage2.py writes all of that straight out of the run folders.
"""

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

# How close the measured ratio has to sit to the closed form. Two percent, the
# same band Stage 1 declared, and for a better reason: there is no model-form
# difference hiding in this comparison, so the band really is a numerical one.
BAND_PCT = 2.0

# Sliding interface energy screening criterion, percent of internal energy.
# It applies to the frictionless runs only. There the sliding interface energy
# is nothing but penalty energy, so it is a straight measure of how hard the
# contact is working against itself. With friction switched on the same counter
# also collects the frictional dissipation, which is real work and is large -
# at 30 deg and mu = 0.3 it comes out near 80 % of the internal energy. Screening
# that at 5 % would fail a run for doing its job. Those cases are checked against
# friction_work() instead.
SLIDING_LIMIT_PCT = 5.0

# How far total energy / initial energy may drift from 1. This one does apply to
# every case: friction is inside the balance, so a contact that creates energy
# shows up here whether mu is zero or not.
ENERGY_BALANCE_PCT = 1.0


def wedge_ratio(angle_deg, mu):
    """Closed form, repeated here so post-processing does not import the model."""
    ta = math.tan(math.radians(angle_deg))
    denom = 1.0 - mu * ta
    if denom <= 1e-12:
        return math.inf
    return (ta + mu) / denom


def friction_work(angle_deg, mu, p_force, y_tip):
    """
    Work the friction should have done over the stroke, N mm.

    The plate sets the deflection by geometry: once contact is made the tip goes
    down by s tan(a) for a travel s, so the contact force is proportional to s
    for a linear lance and the integral over the stroke is half the end value.
    Stage 1 put the geometric non-linearity at about 1 %, so half the end value
    is good to about that.

        N     = P / (cos a - mu sin a)      normal force at the end of the ramp
        s     = y_tip / tan a               travel from first contact to the end
        slip  = s / cos a                   how far the tip slides along the face

        E     = 0.5 mu N slip

    Measured sliding interface energy also carries the penalty energy, which the
    frictionless runs show at around 2 % of internal energy, so the measured
    value should sit a little above this.
    """
    a = math.radians(angle_deg)
    denom = math.cos(a) - mu * math.sin(a)
    if denom <= 0.0:
        return float("nan")
    normal = p_force / denom
    slip = y_tip / math.tan(a) / math.cos(a)
    return 0.5 * mu * normal * slip


def load_targets():
    path = RESULTS / "analytical_targets.json"
    if not path.exists():
        raise SystemExit("Run scripts/analytical.py first - targets are missing.")
    return json.loads(path.read_text())


def measured_ratio(r):
    """
    W/P for one run: the mean axial reaction divided by the mean transverse one,
    both averaged over the sliding window.

    Two separate decisions are buried in that sentence, and both were got wrong
    first time and corrected against data.

    WHICH WINDOW. W/P = (tan a + mu)/(1 - mu tan a) is a sliding relation - it
    comes from Coulomb friction with relative motion at the interface. The first
    version of this script read it during the hold, because that is where the
    force is flat. But during the hold the plate is stationary and the contact
    is stuck, and a stuck penalty interface carries whatever tangential force
    its stick spring held when motion stopped, less whatever it sheds as the
    lance settles. Read there, the four friction runs come out 0.89, 1.52, 1.90
    and 0.57 % low, at roughly four standard errors each, with the miss growing
    with mu - while both frictionless runs sit at 0.03 and -0.07 %, well inside
    their own noise. With mu = 0 there is no tangential force to relax, which is
    exactly why those two are untouched. Read while the plate is moving, all six
    land on the closed form. The frictional work, which is integrated over the
    moving part of the stroke, had already said the same thing from the other
    side: it comes out above the closed form, not below, so friction was fully
    mobilised while sliding.

    The window is not delicate. Over all 35 windows between 0.8 and 1.5 ms long
    inside 2.2 to 4.2 ms, every case has a median deviation within 0.2 % and no
    single window in any case exceeds 0.93 %. 2.5 to 3.5 ms is simply the
    fastest part of the stroke with the force already at half its final value.

    RATIO OF MEANS, not the per-sample ratio averaged over the window. The
    first version of this script used that instead, on the argument that the
    noise is a fluctuation in the contact force magnitude and would therefore
    divide out. The a45 mu=0 run settled the question and the argument was wrong:

      - W scatters 11.7 %, P scatters 5.8 %. A common magnitude fluctuation
        would give both the same relative scatter, so that is not what this is.
        Almost all of it sits in the axial reaction, where the lance is far
        stiffer than in bending and the mode rings well above the 1e-5 s output
        interval - consecutive samples alternate high-low in W and barely at all
        in P.
      - ratio of means gave 0.9993 against a true value of exactly 1. Mean of
        ratios gave 0.9966. Dividing by a fluctuating denominator sample by
        sample biases the result; the ratio of means is the impulse ratio over
        the window and does not.

    Both corrections came from the same run and point the same way: the estimator
    has to be unbiased and it has to be read in the state the relation describes.
    """
    return r["W"] / r["P"]


def load_results():
    path = RESULTS / "stage2_results.csv"
    if not path.exists():
        raise SystemExit(f"Fill in {path.relative_to(REPO)} from the runs first.")

    rows = []
    with open(path) as f:
        for row in csv.DictReader(line for line in f if not line.startswith("#")):
            if not row["W_N"].strip():
                continue                      # not run yet
            rows.append({
                "case": row["case"],
                "angle": float(row["angle_deg"]),
                "mu": float(row["mu"]),
                "disp": abs(float(row.get("tip_disp_mm") or "nan")),
                "W": abs(float(row["W_N"])),
                "P": abs(float(row["P_N"])),
                "W_hold": abs(float(row.get("W_hold_N") or "nan")),
                "P_hold": abs(float(row.get("P_hold_N") or "nan")),
                # Recorded for the report, not used in the comparison.
                "scatter": float(row.get("ratio_scatter_pct") or "nan"),
                "ke_ie": float(row.get("ke_over_ie_pct") or "nan"),
                "sl_ie": float(row.get("sliding_over_ie_pct") or "nan"),
                "ie": float(row.get("internal_energy_Nmm") or "nan"),
                "eratio": float(row.get("energy_ratio") or "nan"),
            })
    if not rows:
        raise SystemExit("No results filled in yet.")
    return rows


def plot_friction(rows, angles):
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    colours = {30.0: "#9d4c19", 45.0: "#2f6f6a"}

    mus = [i / 200.0 for i in range(0, 71)]          # 0 to 0.35
    for a in angles:
        col = colours.get(a, "#5a646d")
        ax.plot(mus, [wedge_ratio(a, m) for m in mus], "-", color=col, lw=1.4,
                label=f"closed form, {a:.0f}°")
        pts = sorted((r for r in rows if r["angle"] == a), key=lambda r: r["mu"])
        if pts:
            ax.plot([p["mu"] for p in pts], [measured_ratio(p) for p in pts],
                    "o", color=col, ms=7, mfc="white", mew=1.8,
                    label=f"LS-DYNA, {a:.0f}°")

    ax.set_xlabel("friction coefficient  µ")
    ax.set_ylabel("W / P   [-]")
    ax.set_title("Stage 2 wedge relation, axial over transverse force")
    ax.set_xlim(-0.01, 0.355)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(fontsize=8, loc="best", framealpha=0.92)
    fig.tight_layout()

    out = RESULTS / "figures" / "stage2_wedge_ratio.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


def plot_history(rows):
    """
    W/P against time for every run, against its own closed form.

    This is the figure the window argument rests on: flat along the closed form
    while the plate moves, sagging once it stops, and sagging further the larger
    mu is - while the two mu = 0 runs stay flat all the way through.

    Plotted as a running 0.5 ms average, and as a ratio of the two averages
    rather than an average of the ratios - the same estimator the table uses,
    for the same reason. The per-sample ratio is unreadable: the lance's axial
    mode rings well above the 1e-5 s output interval and swings it several
    percent between consecutive samples. That is noise on the measurement, not
    behaviour of the contact, and it is what the averaging is there to remove.
    """
    path = RESULTS / "stage2_ratio_history.csv"
    if not path.exists():
        return None

    raw = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            raw.setdefault(row["case"], []).append(
                (float(row["time_ms"]), abs(float(row["W_N"])), abs(float(row["P_N"]))))

    half = 25                                   # 0.5 ms at a 0.01 ms interval
    series = {}
    for case, pts in raw.items():
        pts.sort()
        out = []
        for i in range(half, len(pts) - half):
            block = pts[i - half:i + half + 1]
            p = sum(b[2] for b in block)
            if p <= 0.0:
                continue
            out.append((pts[i][0], sum(b[1] for b in block) / p))
        series[case] = out

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colours = {30.0: "#9d4c19", 45.0: "#2f6f6a"}
    styles = {0.00: ":", 0.10: "-.", 0.20: "--", 0.30: "-"}

    for r in rows:
        pts = series.get(r["case"])
        if not pts:
            continue
        col = colours.get(r["angle"], "#5a646d")
        target = wedge_ratio(r["angle"], r["mu"])
        # Normalised, so six runs with targets from 0.58 to 1.50 share one axis.
        ax.plot([t for t, _ in pts], [v / target for _, v in pts],
                styles.get(r["mu"], "-"), color=col, lw=1.1,
                label=f"{r['angle']:.0f}°, µ = {r['mu']:.2f}")

    ax.axhline(1.0, color="#333333", lw=1.0)
    ax.axvspan(2.5, 3.5, color="#4a7c59", alpha=0.10)
    ax.axvspan(5.0, 5.5, color="#8a8a8a", alpha=0.12)
    ax.text(3.0, 0.9745, "measured here\nplate sliding", ha="center", fontsize=7.5)
    ax.text(5.25, 0.9745, "plate\nstopped", ha="center", fontsize=7.5)
    ax.set_xlabel("time  [ms]")
    ax.set_ylabel("W / P   divided by its closed form")
    ax.set_title("Stage 2, W/P through the stroke, running 0.5 ms average")
    ax.set_xlim(1.5, 5.4)
    ax.set_ylim(0.970, 1.030)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(fontsize=7.5, loc="upper right", ncol=2, framealpha=0.92)
    fig.tight_layout()

    out = RESULTS / "figures" / "stage2_ratio_history.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


def main():
    everything = load_targets()
    y = everything["geometry"]["y"]
    rows = load_results()

    print("Stage 2 - contact and friction, checked on W/P\n")
    print("The ratio carries no stiffness, so the Stage 1 deviation from beam")
    print("theory does not enter this comparison.\n")

    # The deflection is checked first, as in Stage 1, but it is not required to
    # be exact here and it is not a pass/fail. The plate face is inclined, so
    # the height of the face under the tip depends on where the tip corner has
    # moved to, and the corner both draws in and rides forward as the lance
    # bends. A couple of percent either side of 0.600 is the geometry, not an
    # error. It is printed because a tip nowhere near 0.600 would mean the
    # forces were read at the wrong place in the stroke.
    print("Tip deflection reached")
    for r in rows:
        if math.isnan(r["disp"]):
            print(f"  {r['case']:>12}  not recorded")
            continue
        off = (r["disp"] / y - 1.0) * 100.0
        flag = "" if abs(off) < 5.0 else "   <-- nowhere near the stroke end, check the run"
        print(f"  {r['case']:>12}  {r['disp']:.4f} mm  {off:+.2f} %{flag}")

    print(f"\n{'case':>12} {'angle':>6} {'mu':>6} {'W [N]':>9} {'P [N]':>9} "
          f"{'W/P':>8} {'target':>8} {'dev':>8} {'scatter%':>9}")
    worst, worst_case = 0.0, None
    for r in rows:
        ratio = measured_ratio(r)
        target = wedge_ratio(r["angle"], r["mu"])
        dev = (ratio / target - 1.0) * 100.0
        if abs(dev) > abs(worst):
            worst, worst_case = dev, r["case"]
        flag = "" if abs(dev) <= BAND_PCT else "  outside band"
        print(f"{r['case']:>12} {r['angle']:>5.0f}° {r['mu']:>6.2f} "
              f"{r['W']:>9.4f} {r['P']:>9.4f} {ratio:>8.4f} {target:>8.4f} "
              f"{dev:>7.2f}% {r['scatter']:>9.1f}{flag}")

    # The regression test. At the retention angle without friction the relation
    # reduces to W = P, and nothing but the contact normal decides it: no
    # friction model, no stiffness, no material. If this one is wrong, the
    # contact is wrong and every other number in the table is unsafe.
    reg = [r for r in rows if abs(r["angle"] - 45.0) < 1e-9 and r["mu"] == 0.0]
    print("\nRegression test - 45 deg, frictionless, W/P must be 1.000")
    if not reg:
        print("  Not run yet. Run this one first; it is the cheapest way to find")
        print("  out whether the contact formulation is doing what it should.")
    else:
        r = reg[0]
        ratio = measured_ratio(r)
        print(f"  measured {ratio:.4f}   deviation {(ratio - 1.0) * 100:+.2f} %")
        print("  Nothing but the contact normal sets this number - no friction,")
        print("  no stiffness, no material. It is the strongest single check in")
        print("  Stage 2.")

    # Frictionless cases isolate the contact formulation from the friction model.
    free = [r for r in rows if r["mu"] == 0.0]
    if free:
        print("\nContact formulation, frictionless cases only")
        for r in free:
            ratio = measured_ratio(r)
            target = wedge_ratio(r["angle"], 0.0)
            print(f"  {r['angle']:>3.0f}°  {ratio:.4f} vs tan = {target:.4f}"
                  f"   {(ratio / target - 1.0) * 100:+.2f} %")
        print("  With mu = 0 the ratio is pure geometry. A deviation here is the")
        print("  contact algorithm, not the friction coefficient.")

    print("\nSolution quality")
    print(f"{'case':>12} {'KE/IE %':>9} {'sliding/IE %':>13} {'energy ratio':>14}")
    for r in rows:
        flag = ""
        if r["mu"] == 0.0 and r["sl_ie"] > SLIDING_LIMIT_PCT:
            flag = "   <-- penalty energy too high"
        if abs(r["eratio"] - 1.0) * 100.0 > ENERGY_BALANCE_PCT:
            flag = "   <-- energy balance broken"
        print(f"{r['case']:>12} {r['ke_ie']:>9.3f} {r['sl_ie']:>13.3f} "
              f"{r['eratio']:>14.5f}{flag}")
    print(f"  KE/IE is the max over the second half of the ramp. Earlier than")
    print(f"  that the lance holds almost no internal energy and the ratio is")
    print(f"  set by its own denominator.")
    print(f"  Energy ratio is screened at ±{ENERGY_BALANCE_PCT:.0f} % for every case.")
    print(f"  Sliding interface energy is screened at {SLIDING_LIMIT_PCT:.0f} % of")
    print(f"  internal for the frictionless runs only - see below.")

    # With friction on, the sliding interface energy is mostly frictional work,
    # so it is compared with the closed form rather than screened. This is a
    # second, independent check on the same friction law the ratio tests: the
    # ratio is a force statement, this one is an energy statement.
    rubbing = [r for r in rows if r["mu"] > 0.0 and not math.isnan(r["ie"])]
    if rubbing:
        print("\nFrictional work, measured against closed form")
        print(f"{'case':>12} {'measured':>10} {'closed form':>12} {'diff':>8}")
        for r in rubbing:
            meas = r["sl_ie"] / 100.0 * r["ie"]
            # The end-of-stroke transverse force, not the sliding-window one:
            # the work is an integral over the whole stroke and is normalised on
            # the force the stroke finishes at.
            calc = friction_work(r["angle"], r["mu"], r["P_hold"], y)
            d = (meas / calc - 1.0) * 100.0 if calc else float("nan")
            print(f"{r['case']:>12} {meas:>10.4f} {calc:>12.4f} {d:>7.1f}%")
        print("  N mm at the end of the ramp. The closed form is first order: it")
        print("  takes the contact force as linear in the plate travel and the")
        print("  slip from rigid-body geometry. Everything it leaves out adds to")
        print("  the measured side - the penalty energy the frictionless runs put")
        print("  at about 2 % of internal, the extra slip as the tip corner rides")
        print("  forward with the section rotation, and every small reversal at")
        print("  the interface while the lance rings, which dissipates whichever")
        print("  way it goes. So the measured value sitting above the closed form")
        print("  by around 10 % is the expected sign, and it is what says friction")
        print("  was fully mobilised while the plate was moving.")

    # The same ratio read where the plate has stopped. This is the finding of the
    # stage, so it is printed rather than quietly dropped.
    stopped = [r for r in rows if not math.isnan(r["P_hold"])]
    if stopped:
        print("\nThe same ratio read while the contact is stuck, hold 5.0-5.5 ms")
        print(f"{'case':>12} {'mu':>6} {'sliding':>9} {'stuck':>9} {'target':>9} "
              f"{'sliding':>9} {'stuck':>9}")
        for r in stopped:
            target = wedge_ratio(r["angle"], r["mu"])
            slide = measured_ratio(r)
            stuck = r["W_hold"] / r["P_hold"]
            print(f"{r['case']:>12} {r['mu']:>6.2f} {slide:>9.4f} {stuck:>9.4f} "
                  f"{target:>9.4f} {(slide / target - 1) * 100:>8.2f}% "
                  f"{(stuck / target - 1) * 100:>8.2f}%")
        print("  W/P is a sliding relation. Once the plate stops there is no")
        print("  relative motion, the contact sticks, and the tangential force")
        print("  relaxes off the stick spring. The two mu = 0 runs are untouched,")
        print("  because with no friction there is nothing to relax - which is")
        print("  what identifies the mechanism rather than leaving it a guess.")

    print("\nSummary")
    if worst_case:
        print(f"  worst deviation from the closed form   {worst:+.2f} %  ({worst_case})")
        print(f"  band                                   ±{BAND_PCT:.0f} %")
    angles = sorted({r["angle"] for r in rows})
    if len(rows) >= 2:
        print("\nFigures")
        print(f"  {plot_friction(rows, angles).relative_to(REPO)}")
        hist = plot_history(rows)
        if hist:
            print(f"  {hist.relative_to(REPO)}")


if __name__ == "__main__":
    main()
