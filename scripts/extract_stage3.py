"""
Reads the Stage 3 run folders and fills results/stage3_results.csv.

    python scripts/extract_stage3.py <folder holding the runs>

Needs results/stage3_model.json from make_stage3_models.py for the node
coordinates, the tooth stations and each case's stroke and baseline.

What is measured, and why these quantities
------------------------------------------
The pull-out force is the summed AXIAL reaction at the root, not the transverse
one. Lance equilibrium in x settles it: the TPA is frictionless against a flat
horizontal face so it puts no x-force into the lance at all, which leaves the
root's x-reaction equal and opposite to the terminal's axial load. So

    W = |sum Fx at the root|        the pull-out force, Stage 2's W
    P = |sum Fz at the root|        the transverse force, Stage 2's P

and for the runs with no TPA, or before the TPA engages, W/P is exactly the
Stage 2 wedge ratio. That lets the effective face angle be read straight out of
the model by inverting it:

    W/P = (tan b + mu)/(1 - mu tan b)   ->   tan b = (W/P - mu)/(1 + mu W/P)

At mu = 0 it is just b = atan(W/P), with no friction coefficient in it at all.

The same angle is obtained a second, independent way - from the deflected
shape. Two spine nodes either side of the retention face give the lance's local
rotation directly, and the effective angle is then the drawn angle minus that
rotation. The two routes share no data: one is forces, the other is geometry.

Where to read them, and why it is not one instant
-------------------------------------------------
The rotation grows with the deflection, so the effective angle is not a number
but a sweep: it starts at the drawn angle and shallows as the lance bends -
about 42 deg early in the travel down to about 37.3 at release. Any single
value only means something if the instant is named, and comparing two routes at
one instant tests far less than comparing both sweeps.

Two earlier versions got this wrong. The first read the raw peak, which is the
ringing amplitude and falls at whatever travel the ring happens to crest at.
The second read a filtered value at release, a point reading on a signal too
noisy to support one - the answer then moved 7 deg with the filter cutoff.

Forces are filtered, displacements are not
------------------------------------------
The contact force rings far above anything the loading contains. The first set
of runs wrote it every 1e-5 s, too slow to resolve that ring, so it aliased -
and aliased data cannot be filtered clean afterwards. It showed up as an
apparent frequency of 44 kHz in one run and 16 kHz in the same model run at
half the rate, which no real structural mode can do.

So forces are written fast enough to resolve the ring, and samples_per_cycle
reports whether that worked.

Resolving it is necessary and not sufficient. The ripple is 55 to 72 % of the
force's own level, and W and P are two projections of the same noisy contact
impulse, so their ratio at an instant scatters by several degrees however it is
filtered. The reading is therefore the ratio of the MEANS over a window, at
several stations through the stroke, each compared against the rotation over
the same window - which is Stage 2's estimator, applied to a quantity that is
moving. Two columns say whether it held: angle_gap_max_deg, the worst
disagreement between the two routes, and window_spread_deg, which says whether
the answer depends on the window.

Displacements stay at the slower rate and are not filtered. A displacement is
the double integral of the acceleration, so the ring is already small in it -
which is why the rotation measurement came through the first set intact while
every force number in it had to be thrown away.
"""

import csv
import json
import math
import re
import sys
from pathlib import Path

from extract_stage2 import read_glstat, nearest, mean

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
FLOAT = re.compile(r"[-+]?\d*\.\d+[Ee][-+]\d+")

T_RAMP, T_END = 5.0e-3, 5.5e-3
T_RAMP_INSERT, T_END_INSERT = 30.0e-3, 31.0e-3
TOL = 1.0e-6

# Reporting cutoff, and the three others the answer is repeated at. The band
# that matters is set by the loading: a 5 ms quintic ramp puts essentially
# nothing above a few hundred Hz, and the lance's first bending mode is
# 5.1 kHz. A 10 kHz cutoff keeps both and removes the contact ring above them.
# The sweep is what makes the choice defensible - if the answer moves with the
# cutoff, the cutoff is answering.
CUTOFF_HZ = 10.0e3
MIN_SAMPLES_PER_CYCLE = 8.0       # below this the force output is not resolved
MIN_RIPPLE = 0.05                 # ...but only if the ripple is worth resolving

# How the force is read, and why it is not read at an instant.
#
# The contact force ripples by 55 to 72 % of its own level. W and P are two
# projections of the same noisy contact impulse, so their ratio AT A POINT has
# several degrees of scatter in it no matter how the signal is filtered -
# filtering a point reading just trades one arbitrary number for another, and
# the answer then moves with the filter.
#
# Stage 2 settled this and the rule is its rule: take the ratio of the MEANS
# over a window, not the mean of the ratio, and not a point. The window has to
# be long against the ringing and short against the stroke. The ring is above
# 100 kHz (under 10 us) and the ramp is 5 ms, so a few hundred microseconds
# sits cleanly between them, and 150 us also covers most of a cycle of the
# lance's own 5.1 kHz bending mode, which is dynamic overshoot rather than the
# quasi-static force.
#
# The effective angle changes through the stroke as the lance rotates, so a
# single window would average an angle that is moving. The reading is therefore
# taken at several stations and compared against the rotation measured over the
# SAME window at the SAME station. That is the real two-route check: not one
# number against one number, but the whole sweep against the whole sweep.
# Where the stations go, and why the variable is the lance's lift and not the
# terminal's travel.
#
# The wedge relation describes ONE face. The insertion run crosses three, and
# the first attempt at this placed stations by travel fraction and walked into
# two of them:
#
#   lift below the terminal's nose chamfer depth
#           the contact is on the chamfer, which is drawn at 45 deg, so the
#           force route reads about 35.5 deg against a shape route computing
#           30 deg plus the rotation. Off by 4 deg, and correctly so - it is
#           measuring a different surface.
#   lift at the full protrusion
#           the corner has run off the lead-in onto the flat crest. The face
#           angle goes to zero, W/P collapses towards mu, and the force route
#           reads 15.7 deg against 37.3. Off by 21 deg, again correctly.
#   in between
#           the contact is on the lead-in and the two routes agree to 0.4 deg.
#
# Both boundaries are geometry and both are known from the deck before any run:
# the chamfer depth and the tooth protrusion. So the progress variable is the
# crest lift, which is what actually says where on the tooth the contact is,
# and the stations span the middle of it. Travel does not say that - it is the
# same number whichever surface is carrying.
STATIONS = (0.35, 0.50, 0.65, 0.80, 0.90)   # fractions of the usable lift span
WINDOW_US = 150.0
WINDOWS_US = (50.0, 100.0, 150.0, 250.0, 400.0)
GAP_LIMIT = 1.5                   # degrees, force against shape, worst station
WINDOW_SPREAD_LIMIT = 1.5         # degrees across the window sweep


def effective_from_shape(case, rotation_deg):
    """
    The face angle the terminal actually meets, from the drawn angle and the
    lance's measured rotation. Both angles are from the insertion axis.

    The sign is opposite on the two faces and that is the whole point of the
    stage. The lead-in runs uphill towards the tip and the retention face runs
    downhill towards the root, so one nose-up rotation tilts them in opposite
    senses: the lead-in gets STEEPER by theta, the retention gets SHALLOWER by
    the same theta. A 30 deg face and a 45 deg face converge on each other, and
    the asymmetry the pair was drawn to produce largely goes away.

    Getting this sign wrong makes the insertion run disagree with its own
    prediction by twice the rotation while still looking self-consistent,
    because the force route is derived from the same geometry.
    """
    if case.startswith("insert"):
        return 30.0 + abs(rotation_deg)
    return 45.0 - abs(rotation_deg)


def smoothstep(s):
    s = min(max(s, 0.0), 1.0)
    return 10 * s**3 - 15 * s**4 + 6 * s**5


def read_spcforc_nodes(path, z_of):
    """[(time, sum Fx, sum Fz, sum z*Fx)] - one entry per output block."""
    out, t, fx, fz, m, started = [], None, 0.0, 0.0, 0.0, False
    with open(path, errors="ignore") as f:
        for line in f:
            if "output at time" in line:
                if started:
                    out.append((t, fx, fz, m))
                mt = FLOAT.search(line)
                t, fx, fz, m, started = (float(mt.group()) if mt else None,
                                         0.0, 0.0, 0.0, True)
            elif line.lstrip().startswith("node="):
                nid = int(re.match(r"\s*node=\s*(\d+)", line).group(1))
                v = FLOAT.findall(line)
                if len(v) < 3 or nid not in z_of:
                    continue
                fx += float(v[0])
                fz += float(v[2])
                m += z_of[nid] * float(v[0])
    if started:
        out.append((t, fx, fz, m))
    return out


def read_nodout_nodes(path, wanted):
    """{node id: [(time, z displacement, current x)]}."""
    out = {n: [] for n in wanted}
    t, in_block = None, False
    with open(path, errors="ignore") as f:
        for line in f:
            if "n o d a l" in line:
                mt = re.search(r"at time\s+([-+]?\d*\.\d+[Ee][-+]\d+)", line)
                t, in_block = (float(mt.group(1)) if mt else None), False
            elif "nodal point" in line:
                in_block = True
            elif in_block:
                mt = re.match(r"\s*(\d+)\s", line)
                if not mt:
                    in_block = False
                    continue
                nid = int(mt.group(1))
                if nid in out:
                    v = [float(x) for x in FLOAT.findall(line)]
                    if len(v) >= 12:
                        out[nid].append((t, v[2], v[9]))
    return out


def find_run(root, case, all_cases=()):
    """
    The folder holding one case's output.

    Exact name first. Stage 2's version matched on substring alone, which is
    fine there and wrong here: extract_mu020 is contained in
    extract_mu020_fine and extract_mu020_slow, so a plain substring search can
    hand the baseline the sensitivity run's results and then compare that run
    against itself.

    The substring fallback is kept, for a folder named by hand or with a
    suffix, but it refuses two things: a folder that also matches some other
    case, and an ambiguous match. Returning nothing is recoverable. Returning
    the wrong run is not - it produces a plausible number with no sign that
    anything went wrong.
    """
    if not root.is_dir():
        return None
    exact = {f"stage3_{case}", case}
    candidates = [root] + [d for d in sorted(root.rglob("*")) if d.is_dir()]
    runs = [d for d in candidates if (d / "spcforc").exists()]
    for d in runs:
        if d.name in exact:
            return d

    others = [c for c in all_cases if c != case and case in c]
    loose = [d for d in runs
             if case in d.name and not any(o in d.name for o in others)]
    return loose[0] if len(loose) == 1 else None


def node_maps(model, case):
    """
    The node-to-coordinate maps for one case.

    A case that builds its own mesh carries its own map. Using the shared one
    for it reads every displacement at the wrong place, and does it silently,
    because the root nodes are numbered first and keep their IDs - so the
    forces stay right while the geometry quietly does not.
    """
    m = model["cases"][case].get("nodes") or model
    return ({int(k): v for k, v in m["root_nodes"].items()},
            {int(k): v for k, v in m["spine_nodes"].items()},
            {int(k): v for k, v in m["under_nodes"].items()})


def cfc(series, dt, cutoff):
    """
    Zero-phase 2-pole Butterworth low pass, the SAE J211 filter.

    Applied forward then backward, so it has no phase shift and cannot move
    the peak in time. Written out rather than imported: this runs on whatever
    Python is next to LS-DYNA, and a filter is not worth a dependency.
    """
    wd = 2.0 * math.pi * cutoff * 2.0775
    wa = math.tan(wd * dt / 2.0)
    den = 1.0 + math.sqrt(2.0) * wa + wa * wa
    a0 = wa * wa / den
    a1, a2 = 2.0 * a0, a0
    b1 = -2.0 * (wa * wa - 1.0) / den
    b2 = (-1.0 + math.sqrt(2.0) * wa - wa * wa) / den

    def pass_once(x):
        y = list(x[:2])
        for i in range(2, len(x)):
            y.append(a0 * x[i] + a1 * x[i - 1] + a2 * x[i - 2]
                     + b1 * y[i - 1] + b2 * y[i - 2])
        return y

    if len(series) < 5:
        return list(series)
    return pass_once(pass_once(series)[::-1])[::-1]


def samples_per_cycle(series, dt):
    """
    How many samples the output writes per cycle of the fastest thing in it.

    Returns (samples per cycle, apparent frequency in Hz). The estimate is the
    lag-1 autocorrelation of the signal about its slow trend: an oscillation
    resolved by many samples has neighbouring samples almost equal, r -> +1,
    while one at the Nyquist limit alternates sign every sample, r -> -1. The
    cycle length follows as 2*pi/acos(r).

    This check exists because the first Stage 3 set was written every 1e-5 s
    and every force number in it turned out to be aliased. Nothing in the
    output said so - the curves looked like noisy curves. What gave it away
    was that the apparent frequency differed between two runs of the same
    model, 44 kHz against 16 kHz, which a real structural mode cannot do.

    Below about 8 samples per cycle the signal is not resolved, filtering it
    is meaningless, and the answer will depend on the filter rather than on
    the model. The threshold is stated, not tuned: it is where a Butterworth
    low pass still has a usable transition band under the Nyquist limit.

    Also returns the ripple as a fraction of the signal level, because being
    unresolved only matters if there is something there. Any run has a thin
    white floor from round-off, and white noise is unresolved by definition -
    it alternates every sample. Gating on the cycle count alone would fail a
    clean run for having a 0.3 % wobble on it, which is how the regression run
    first tripped this check.
    """
    n = len(series)
    if n < 32 or dt <= 0:
        return float("nan"), float("nan"), float("nan")
    # Take the trend out first, or a small ripple riding on a large ramp is
    # hidden by the ramp - which is the case here, since the contact force
    # climbs through the whole stroke.
    w = max(8, min(65, n // 8)) | 1
    h = w // 2
    d = []
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        d.append(series[i] - sum(series[a:b]) / (b - a))
    v = sum(x * x for x in d)
    if v <= 0:
        return float("inf"), 0.0, 0.0
    r = sum(a * b for a, b in zip(d, d[1:])) / v
    r = min(max(r, -0.999999), 0.999999)
    spc = 2.0 * math.pi / math.acos(r)
    level = sum(abs(x) for x in series) / n
    ripple = math.sqrt(v / n) / level if level > 0 else 0.0
    return spc, 1.0 / (spc * dt), ripple


def effective_angle(ratio, mu):
    """
    Invert the wedge relation for the face angle. [degrees]

        tan b = (W/P - mu) / (1 + mu W/P)

    The friction coefficient is an input to the model, so this is not circular:
    the only measured quantity is the force ratio, and at mu = 0 not even that
    - the angle is then the arctangent of the ratio and nothing else.
    """
    if ratio <= 0.0 or not math.isfinite(ratio):
        return float("nan")
    denom = 1.0 + mu * ratio
    return math.degrees(math.atan((ratio - mu) / denom)) if denom > 0 else float("nan")


def extract_lance_only(folder, model, case, spec):
    """
    The regression run: the tip is lifted by a prescribed amount and nothing
    else touches the lance, so the root's transverse reaction is the Stage 1
    tip force with the tooth's stiffening added. If this one drifts, nothing
    downstream of it means anything, because the tooth has changed the beam
    into something Stage 1 no longer describes.
    """
    row, notes = {"case": case, "mu": spec["mu"],
                  "t_ramp_ms": spec["t_ramp"] * 1e3,
                  "tooth_h_mm": spec.get("tooth_h")}, []
    t_ramp = spec["t_ramp"]

    z_of = {int(k): v for k, v in model["root_nodes"].items()}
    under = {int(k): v for k, v in model["under_nodes"].items()}
    n_tip = min(under, key=lambda n: abs(under[n] - model["tooth_stations"][3]))

    forces = read_spcforc_nodes(folder / "spcforc", z_of)
    traces = read_nodout_nodes(folder / "nodout", {n_tip})
    if not forces:
        notes.append("spcforc empty")
        return row, notes

    t = [f[0] for f in forces]
    dt = (t[-1] - t[0]) / max(1, len(t) - 1)
    row["force_dt_us"] = dt * 1e6
    band = [i for i, ti in enumerate(t) if 0.2 * t_ramp <= ti <= t_ramp]
    if len(band) > 64:
        spc, f_app, ripple = samples_per_cycle([abs(forces[i][2]) for i in band], dt)
        row["samples_per_cycle"] = spc
        row["ringing_kHz"] = f_app / 1000.0
        row["ripple_pct"] = ripple * 100.0

    # No contact in this run, so the reaction is smooth and the filter has
    # almost nothing to remove. It is applied anyway, at the same cutoff as
    # everywhere else, so the regression is measured the way the rest is.
    row["cutoff_Hz"] = CUTOFF_HZ
    Pf = cfc([abs(f[2]) for f in forces], dt, CUTOFF_HZ)
    i = min(range(len(t)), key=lambda k: abs(t[k] - t_ramp))
    if abs(t[i] - t_ramp) > 0.05 * t_ramp:
        notes.append(f"no output near the end of the ramp - closest is "
                     f"{t[i] * 1e3:.3f} ms against {t_ramp * 1e3:.3f} ms")
    row["P_N"] = Pf[i]
    row["W_N"] = abs(forces[i][1])
    row["W_raw_peak_N"] = max(abs(f[2]) for f in forces)
    row["read_at"] = "at the end of the ramp"
    row["read_travel_mm"] = 0.0

    lift = None
    if traces.get(n_tip):
        lift = min(traces[n_tip], key=lambda s: abs(s[0] - t_ramp))[1]
        row["tip_lift_mm"] = abs(lift)

    pred = spec.get("P_predicted")
    if pred:
        dev = 100.0 * (row["P_N"] / pred - 1.0)
        notes.append(f"root reaction {row['P_N']:.4f} N against "
                     f"{pred:.4f} N predicted, {dev:+.2f} %")
        notes.append(f"baseline: {spec.get('baseline_is', 'unstated')}")
    if lift is not None and spec.get("tip_lift"):
        notes.append(f"prescribed lift {spec['tip_lift']:.3f} mm, tip moved "
                     f"{abs(lift):.4f} mm")

    row.update(read_quality(folder, t_ramp))
    return row, notes


def read_quality(folder, t_ramp, t_peak=None):
    """
    The glstat columns every run is screened on.

    Read at the peak force, not at the end of the ramp. The two are not the
    same instant here: the lance releases at about 45 % of the stroke and the
    internal energy collapses afterwards, so a ratio taken at the end of the
    ramp has a near-zero denominator and reports the denominator rather than
    the dynamics. That is the same mistake Stage 2 made at the other end of
    the run, and the rule it produced is the one applied here - measure the
    solution quality at the instant the answer is taken from.

    The energy balance is the exception. It is a global conservation check
    over the whole run, so it is read from the last block.
    """
    out = {}
    if not (folder / "glstat").exists():
        return out
    g = read_glstat(folder / "glstat")
    if not g:
        return out
    ref = t_peak if t_peak else t_ramp

    ie_max = max((b.get("internal energy", 0.0) for b in g), default=0.0)
    window = [b for b in g if 0.5 * ref <= b["time"] <= ref + TOL
              and b.get("internal energy", 0.0) > 0.10 * ie_max]
    ke = [b["kinetic energy"] / b["internal energy"] * 100.0 for b in window]
    if ke:
        out["ke_over_ie_pct"] = max(ke)

    at = nearest(g, ref, key=lambda b: b["time"], window=5 * ref / 100)
    if at and at.get("internal energy", 0.0) > 0.0:
        out["internal_energy_Nmm"] = at["internal energy"]
        out["hourglass_over_ie_pct"] = (at.get("hourglass energy", 0.0)
                                        / at["internal energy"] * 100.0)
        out["sliding_over_ie_pct"] = (at.get("sliding interface energy", 0.0)
                                      / at["internal energy"] * 100.0)
    out["energy_ratio"] = g[-1].get("total energy / initial energy", float("nan"))
    return out


def interpolate(series, times):
    """
    A slow-rate displacement trace put onto the fast force time base.

    Forces and displacements are written at different rates on purpose - the
    force needs 1 MHz to be resolved, the displacement does not need anything
    like it. Linear interpolation is right for the displacement because it is
    smooth at this scale; doing it the other way round, decimating the force
    onto the displacement rate, is what produced the aliased first set.
    """
    if not series:
        return [float("nan")] * len(times)
    ts = [s[0] for s in series]
    out, j = [], 0
    for t in times:
        while j + 2 < len(ts) and ts[j + 1] < t:
            j += 1
        t0, t1 = ts[j], ts[min(j + 1, len(ts) - 1)]
        v0, v1 = series[j][1], series[min(j + 1, len(series) - 1)][1]
        out.append(v0 if t1 <= t0 else v0 + (v1 - v0) * (t - t0) / (t1 - t0))
    return out


def extract(folder, model, case, spec):
    row, notes = {"case": case}, []
    mu = spec["mu"]
    gap = spec["tpa_clearance"]
    stroke = spec["stroke"]
    direction = spec["direction"]
    insert = direction < 0
    # The ramp comes from the model index, not from a constant here: the
    # sensitivity run doubles it, and reading it from the deck's own record is
    # the only way the travel axis stays right for both.
    t_ramp = spec.get("t_ramp") or (T_RAMP_INSERT if insert else T_RAMP)

    row.update({"mu": mu, "tpa_gap_mm": gap, "stroke_mm": stroke,
                "t_ramp_ms": t_ramp * 1e3, "tooth_h_mm": spec.get("tooth_h")})

    z_of, spine, under = node_maps(model, case)
    x1, x2, x3, x4 = model["tooth_stations"]

    # Nodes we need: the two spine stations either side of the retention face
    # for the rotation, the crest underside for the release condition, and the
    # tip underside because that is what the TPA stops.
    def nearest_node(table, target):
        return min(table, key=lambda n: abs(table[n] - target))

    n_rot_a = nearest_node(spine, x1)
    n_rot_b = nearest_node(spine, x2)
    n_crest = nearest_node(under, 0.5 * (x2 + x3))
    n_tip = nearest_node(under, x4)

    forces = read_spcforc_nodes(folder / "spcforc", z_of)
    traces = read_nodout_nodes(folder / "nodout", {n_rot_a, n_rot_b, n_crest, n_tip})
    if not forces or not traces.get(n_crest):
        notes.append("spcforc or nodout empty")
        return row, notes

    t = [f[0] for f in forces]
    W_raw = [abs(f[1]) for f in forces]
    P_raw = [abs(f[2]) for f in forces]
    dt = (t[-1] - t[0]) / max(1, len(t) - 1)
    row["force_dt_us"] = dt * 1e6

    # Is the force output fast enough to mean anything? Read over the loaded
    # part of the ramp, where the contact is carrying.
    band = [i for i, ti in enumerate(t) if 0.2 * t_ramp <= ti <= t_ramp]
    if len(band) > 64:
        spc, f_app, ripple = samples_per_cycle([W_raw[i] for i in band], dt)
        row["samples_per_cycle"] = spc
        row["ringing_kHz"] = f_app / 1000.0
        row["ripple_pct"] = ripple * 100.0
        if spc < MIN_SAMPLES_PER_CYCLE and ripple > MIN_RIPPLE:
            notes.append(f"FORCE NOT RESOLVED - {spc:.1f} samples per cycle at "
                         f"{f_app / 1e3:.0f} kHz apparent. The oscillation is at "
                         f"or past the Nyquist limit of this output rate, so it "
                         f"is folded into the data and no filter can take it "
                         f"back out. Force columns below are not usable.")
        else:
            notes.append(f"force resolved: {spc:.1f} samples per cycle, ripple "
                         f"{ripple * 100:.1f} % of the level")

    # Displacements onto the force time base.
    crest = interpolate(traces[n_crest], t)
    tip = interpolate(traces[n_tip], t)
    ra = interpolate(traces[n_rot_a], t)
    rb = interpolate(traces[n_rot_b], t)
    dx_rot = abs(spine[n_rot_b] - spine[n_rot_a])
    theta = [math.degrees(math.atan((b - a) / dx_rot)) for a, b in zip(ra, rb)]
    travel = [abs(stroke) * smoothstep(ti / t_ramp) for ti in t]

    # One filtered pass to take the ring out, then the windowed ratio below.
    # The filter alone is not the measurement - it cannot be, on a signal whose
    # ripple is most of its own level.
    W = cfc(W_raw, dt, CUTOFF_HZ)
    P = cfc(P_raw, dt, CUTOFF_HZ)
    row["cutoff_Hz"] = CUTOFF_HZ

    ramp_i = [i for i, ti in enumerate(t) if ti <= t_ramp + TOL]
    if not ramp_i:
        notes.append("no output inside the ramp")
        return row, notes

    row["_curve"] = [{"t": t[i], "travel": travel[i], "W": W[i], "P": P[i],
                      "W_raw": W_raw[i], "crest": crest[i], "tip": tip[i],
                      "theta": theta[i]}
                     for i in range(0, len(t), max(1, len(t) // 1200))]

    # Release: the crest has risen by the full protrusion. Geometry only, no
    # force in it, and the one event in the run with a sharp definition.
    rel = [i for i in range(len(t)) if abs(crest[i]) >= model["protrusion"] - 1e-4]
    row["released"] = "yes" if rel else "no"
    if rel:
        i = rel[0]
        row["release_travel_mm"] = travel[i]
        row["rotation_at_release_deg"] = abs(theta[i])
        row["angle_from_shape_deg"] = effective_from_shape(case, theta[i])
        notes.append(f"released at {travel[i]:.4f} mm, lance rotated "
                     f"{abs(theta[i]):.2f} deg -> effective face "
                     f"{row['angle_from_shape_deg']:.2f} deg")
    else:
        notes.append(f"never released, crest reached "
                     f"{max(abs(c) for c in crest):.4f} of "
                     f"{model['protrusion']:.3f} mm")

    def window(centre_i, half_us):
        lo, hi = t[centre_i] - half_us * 1e-6, t[centre_i] + half_us * 1e-6
        return [i for i in ramp_i if lo <= t[i] <= hi]

    def angle_over(idx):
        """Ratio of the MEANS over the window, inverted for the face angle."""
        if len(idx) < 5:
            return None, None, None
        mW = sum(W[i] for i in idx) / len(idx)
        mP = sum(P[i] for i in idx) / len(idx)
        rot = sum(theta[i] for i in idx) / len(idx)
        if mP <= 1e-6:
            return None, None, None
        r = mW / mP
        return r, effective_angle(r, mu), effective_from_shape(case, rot)

    # The usable window: from the start to release, or to the moment the TPA
    # takes over, whichever comes first. Past either the root reaction is no
    # longer just the terminal's load on one face.
    last = len(ramp_i) - 1
    if rel:
        last = min(last, max(k for k, i in enumerate(ramp_i) if i <= rel[0]))
    if gap is not None:
        hit = [i for i in ramp_i if abs(tip[i]) >= gap - 1e-5]
        if hit:
            last = min(last, max(k for k, i in enumerate(ramp_i) if i <= hit[0]))
    lift_span = abs(crest[ramp_i[last]])
    if lift_span <= 1e-6:
        notes.append("the lance never lifted; no stations")
        lift_span = None

    # The chamfer check. On insertion the terminal's nose chamfer carries the
    # contact until the lance has lifted past its depth, and the chamfer is a
    # different angle from the lead-in, so a station below it measures the
    # wrong surface. The bound is geometry, so it is checked rather than
    # assumed - and it is checked against the deck, not against the answer.
    chamfer = model.get("chamfer_depth")
    if insert and chamfer and lift_span:
        if STATIONS[0] * lift_span < chamfer:
            notes.append(f"first station at {STATIONS[0] * lift_span:.3f} mm of "
                         f"lift is inside the {chamfer:.3f} mm nose chamfer - "
                         f"that station measures the chamfer, not the lead-in")

    stations, gaps = [], []
    if lift_span:
        for frac in STATIONS:
            target = frac * lift_span
            here = [i for i in ramp_i[:last + 1] if abs(crest[i]) >= target]
            if not here:
                continue
            i = here[0]
            r, from_force, from_shape = angle_over(window(i, WINDOW_US))
            if from_force is None or from_force != from_force:
                continue
            stations.append({"case": case, "frac": frac, "travel": travel[i],
                             "lift": abs(crest[i]),
                             "W_over_P": r, "from_force": from_force,
                             "from_shape": from_shape,
                             "gap": from_force - from_shape})
            gaps.append(abs(from_force - from_shape))
    row["_stations"] = stations
    if gaps:
        row["angle_gap_max_deg"] = max(gaps)
        row["angle_gap_mean_deg"] = sum(gaps) / len(gaps)
        last = stations[-1]
        row["W_over_P"] = last["W_over_P"]
        row["angle_from_force_deg"] = last["from_force"]
        row["read_travel_mm"] = last["travel"]
        row["read_at"] = f"windowed, {WINDOW_US:.0f} us, {len(stations)} stations"
        notes.append(f"force against shape over {len(stations)} stations: worst "
                     f"{max(gaps):.2f} deg, mean {sum(gaps) / len(gaps):.2f} deg")

    # Window sensitivity, at the middle station. If the answer moves with the
    # window the window is answering, and the same objection applies as to the
    # filter - so it is measured rather than asserted.
    mid = stations[len(stations) // 2]["travel"] if stations else travel[ramp_i[last]]
    i = min(ramp_i, key=lambda k: abs(travel[k] - mid))
    sweep = [angle_over(window(i, h))[1] for h in WINDOWS_US]
    good = [a for a in sweep if a is not None and a == a]
    if len(good) > 1:
        row["window_spread_deg"] = max(good) - min(good)
        if row["window_spread_deg"] > WINDOW_SPREAD_LIMIT:
            notes.append(f"angle moves {row['window_spread_deg']:.1f} deg across "
                         f"windows of {WINDOWS_US[0]:.0f}-{WINDOWS_US[-1]:.0f} us. "
                         f"The window is answering, not the model.")

    # Force level at the last station, for the retention/insertion comparison.
    if stations:
        idx = window(min(ramp_i, key=lambda k: abs(travel[k] - stations[-1]["travel"])),
                     WINDOW_US)
        row["W_N"] = sum(W[i] for i in idx) / len(idx)
        row["P_N"] = sum(P[i] for i in idx) / len(idx)
        row["crest_lift_mm"] = abs(crest[idx[len(idx) // 2]])
        row["tip_lift_mm"] = abs(tip[idx[len(idx) // 2]])
        row["rotation_at_read_deg"] = abs(theta[idx[len(idx) // 2]])
    row["W_raw_peak_N"] = max(W_raw[i] for i in ramp_i)
    if row.get("W_N", 0) > 1e-9:
        row["ring_over_mean"] = row["W_raw_peak_N"] / row["W_N"]

    # Insertion has no release peak worth quoting - the raw maximum is the
    # terminal's first touch, an impact, and not the insertion force. What
    # matters is the plateau while the lance is riding the lead-in.
    if insert:
        seg = [i for i in ramp_i if 0.15 * abs(stroke) <= travel[i] <= 0.45 * abs(stroke)]
        if seg:
            row["W_plateau_N"] = sum(W[i] for i in seg) / len(seg)
            notes.append(f"insertion plateau {row['W_plateau_N']:.3f} N over "
                         f"{travel[seg[0]]:.2f}-{travel[seg[-1]]:.2f} mm; the raw "
                         f"maximum {row['W_raw_peak_N']:.1f} N is first touch")

    # Where the TPA took over, from the tip meeting the clearance.
    if gap is not None:
        hit = [i for i in range(len(t)) if abs(tip[i]) >= gap - 1e-5]
        if hit:
            row["tpa_engage_travel_mm"] = travel[hit[0]]
            row["W_at_tpa_engage_N"] = W[hit[0]]
            notes.append(f"TPA engaged at {travel[hit[0]]:.4f} mm, W = "
                         f"{W[hit[0]]:.3f} N")
        else:
            notes.append("TPA never touched")

    # Solution quality read where the answer is taken from, not at the end
    # of the ramp - the lance releases at about 45 % of the stroke and the
    # internal energy falls away after that.
    q = min(ramp_i, key=lambda k: abs(travel[k] - mid))
    row.update(read_quality(folder, t_ramp, t[q]))
    return row, notes


HEADER = """\
# Stage 3 results, one row per run. Written by scripts/extract_stage3.py - do
# not edit by hand, rerun the script instead.
#
# The gate thresholds that were in force when this file was written are on the
# limits line below, and the post-processor reads them from there rather than
# importing them. That way the gates reported are the ones actually applied,
# even if the script has been edited since, and the two scripts do not share
# module state - which is how a stale import in a notebook kernel managed to
# pair a new post-processor with an old extractor.
# limits: min_spc={min_spc} min_ripple_pct={min_ripple} gap_limit={gap} \
window_spread_limit={wspread} windows_us={windows}
#
# Forces are low-pass filtered. The raw contact force rings far above the
# frequencies the loading contains, and the first set of runs wrote it too
# slowly to resolve, so the whole set had to be discarded. Two columns say
# whether that is fixed here:
#
# samples_per_cycle     output samples per cycle of the ringing. Below 8 the
#                       signal is aliased and every force column is void.
# angle_gap_max_deg     worst disagreement between the force route and the shape
#                       route, over the stations in stage3_stations.csv. The two
#                       share no data, so this is the real error bar.
# window_spread_deg     how far the measured angle moves across averaging windows
#                       of 50 to 400 us. If it moves, the window is answering.
#
# W_N                   summed AXIAL root reaction, filtered, read at release.
#                       This is the pull-out force: the TPA is frictionless on
#                       a flat face so it puts no x-force into the lance, which
#                       leaves the root x-reaction equal to the terminal load.
# P_N                   summed transverse root reaction at the same instant.
# W_raw_peak_N          the unfiltered maximum, kept so the size of the ring
#                       is on record rather than hidden by the filter.
# W_over_P              the Stage 2 wedge ratio, measured.
# angle_from_force_deg  that ratio inverted for the face angle. At mu = 0 it is
#                       atan(W/P) and contains no assumption at all.
# angle_from_shape_deg  the same angle from the deflected shape - drawn angle
#                       minus the lance rotation at release, measured from two
#                       spine nodes. Shares no data with the column above and
#                       does not depend on the filter.
# released              did the crest rise by the full protrusion.
# W_plateau_N           insertion only. The force while the lance rides the
#                       lead-in. The raw maximum on that run is first touch.
"""

COLUMNS = ["case", "mu", "tpa_gap_mm", "stroke_mm", "tooth_h_mm", "t_ramp_ms",
           "force_dt_us", "samples_per_cycle", "ringing_kHz", "ripple_pct",
           "cutoff_Hz",
           "W_N", "P_N", "W_raw_peak_N", "ring_over_mean",
           "W_over_P", "angle_from_force_deg",
           "angle_gap_max_deg", "angle_gap_mean_deg", "window_spread_deg",
           "angle_from_shape_deg", "rotation_at_release_deg",
           "rotation_at_read_deg", "read_at", "read_travel_mm",
           "crest_lift_mm", "tip_lift_mm",
           "released", "release_travel_mm", "W_plateau_N",
           "tpa_engage_travel_mm", "W_at_tpa_engage_N",
           "ke_over_ie_pct", "hourglass_over_ie_pct",
           "sliding_over_ie_pct", "internal_energy_Nmm", "energy_ratio"]

FMT = {"mu": "{:.2f}", "tpa_gap_mm": "{:.3f}", "stroke_mm": "{:+.3f}",
       "tooth_h_mm": "{:.3f}", "t_ramp_ms": "{:.1f}",
       "force_dt_us": "{:.3f}", "samples_per_cycle": "{:.1f}",
       "ringing_kHz": "{:.1f}", "ripple_pct": "{:.2f}", "cutoff_Hz": "{:.0f}",
       "W_N": "{:.4f}", "P_N": "{:.4f}", "W_raw_peak_N": "{:.4f}",
       "ring_over_mean": "{:.2f}", "W_over_P": "{:.4f}",
       "angle_from_force_deg": "{:.2f}", "angle_gap_max_deg": "{:.2f}",
       "angle_gap_mean_deg": "{:.2f}", "window_spread_deg": "{:.2f}",
       "angle_from_shape_deg": "{:.2f}", "rotation_at_release_deg": "{:.2f}",
       "rotation_at_read_deg": "{:.2f}", "read_travel_mm": "{:.4f}",
       "crest_lift_mm": "{:.4f}", "tip_lift_mm": "{:.4f}",
       "release_travel_mm": "{:.4f}", "W_plateau_N": "{:.4f}",
       "tpa_engage_travel_mm": "{:.4f}", "W_at_tpa_engage_N": "{:.4f}",
       "ke_over_ie_pct": "{:.4f}", "hourglass_over_ie_pct": "{:.4f}",
       "sliding_over_ie_pct": "{:.2f}", "internal_energy_Nmm": "{:.4f}",
       "energy_ratio": "{:.5f}"}


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "ls-dyna" / "stage3"
    if not root.is_dir():
        raise SystemExit(f"{root} is not a folder.")
    index = RESULTS / "stage3_model.json"
    if not index.exists():
        raise SystemExit("Run scripts/make_stage3_models.py first.")
    model = json.loads(index.read_text())

    print(f"Reading runs from {root}\n")
    rows = []
    missing = []
    for case, spec in model["cases"].items():
        folder = find_run(root, case, model["cases"])
        if folder is None:
            print(f"{case:>18}  not found")
            rows.append({"case": case})
            missing.append(case)
            continue
        if spec["direction"] == 0:
            row, notes = extract_lance_only(folder, model, case, spec)
        else:
            row, notes = extract(folder, model, case, spec)
        print(f"{case:>18}  {folder.name}")
        for n in notes:
            print(f"                    {n}")
        if "angle_from_shape_deg" in row:
            print(f"                    effective face at release "
                  f"{row['angle_from_shape_deg']:.2f} deg from the shape")
        print()
        rows.append(row)

    hpath = RESULTS / "stage3_curve_history.csv"
    with open(hpath, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "time_ms", "travel_mm", "W_N", "P_N", "W_raw_N",
                    "crest_lift_mm", "tip_lift_mm", "rotation_deg"])
        for r in rows:
            for c in r.get("_curve", []):
                w.writerow([r["case"], f"{c['t'] * 1e3:.5f}", f"{c['travel']:.6f}",
                            f"{c['W']:.6f}", f"{c['P']:.6f}", f"{c['W_raw']:.6f}",
                            f"{abs(c['crest']):.6f}", f"{abs(c['tip']):.6f}",
                            f"{c['theta']:.5f}"])
    spath = RESULTS / "stage3_stations.csv"
    with open(spath, "w", newline="") as f:
        f.write("# The force route against the shape route, station by station.\n"
                "# Stations are placed by the lance's LIFT, not the terminal's\n"
                "# travel, because the lift is what says which surface is in\n"
                "# contact. On insertion the contact crosses three of them - the\n"
                "# terminal's nose chamfer, the 30 deg lead-in, the flat crest -\n"
                "# and the wedge relation describes only one.\n"
                "#\n"
                "# W_over_P is the ratio of the MEAN axial to the MEAN transverse\n"
                "# root reaction over a window, not a point reading. from_shape is\n"
                "# the drawn angle and the lance rotation over the same window.\n"
                "# The two share no data. gap_deg is the whole error bar.\n")
        w = csv.writer(f)
        w.writerow(["case", "frac_of_lift", "lift_mm", "travel_mm", "W_over_P",
                    "angle_from_force_deg", "angle_from_shape_deg", "gap_deg"])
        for r in rows:
            for st in r.get("_stations", []):
                w.writerow([st["case"], f"{st['frac']:.2f}", f"{st['lift']:.4f}",
                            f"{st['travel']:.4f}", f"{st['W_over_P']:.4f}",
                            f"{st['from_force']:.2f}", f"{st['from_shape']:.2f}",
                            f"{st['gap']:+.2f}"])
    for r in rows:
        r.pop("_curve", None)
        r.pop("_stations", None)

    path = RESULTS / "stage3_results.csv"
    with open(path, "w", newline="") as f:
        f.write(HEADER.format(
            min_spc=MIN_SAMPLES_PER_CYCLE, min_ripple=MIN_RIPPLE * 100,
            gap=GAP_LIMIT, wspread=WINDOW_SPREAD_LIMIT,
            windows=",".join(f"{w:.0f}" for w in WINDOWS_US)))
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: (FMT[c].format(r[c])
                            if c in FMT and r.get(c) is not None
                            and isinstance(r.get(c), (int, float))
                            else ("" if r.get(c) is None else r.get(c, "")))
                        for c in COLUMNS})
    print(f"Written  {path}")
    print(f"Written  {hpath}")
    print(f"Written  {spath}")
    done = sum(1 for r in rows if "P_N" in r)
    print(f"{done} of {len(rows)} runs extracted.")
    if missing:
        print("No spcforc found for: " + ", ".join(missing))
        print(f"Expected a folder named stage3_<case> under {root}.")
    print("Next: python scripts/postprocess_stage3.py")


if __name__ == "__main__":
    main()
