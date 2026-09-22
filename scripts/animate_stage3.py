#!/usr/bin/env python3
"""
The Stage 3 finding as an animation: the lance lifts, the tooth rotates, and the
45 deg retention face and the 30 deg lead-in end up at nearly the same angle.

Everything comes from files already in results/, so nothing has to be rerun:

  stage3_curve_history.csv   lift, rotation and travel of extract_mu020
                             (pull-out) and insert_mu020 (push-in)
  stage3_model.json          tooth stations and protrusion

The terminal outline comes from terminal_geometry() in make_stage3_models.py,
the same function the decks are written from, and the prediction from
effective_angles() in analytical.py.

The face angle is the shape route from extract_stage3.py: 45 minus the rotation
for the retention face, 30 plus the rotation for the lead-in. Each face is shown
from its own run. The two runs lift at different speeds, so the frames are
matched by lift, not by time. The animation stops at release, when the crest has
lifted by the full protrusion.

The lance is drawn to scale, but it is not the FE mesh. From the root to the
tooth it is a cubic through the measured lift and rotation, and the tooth is
moved as one rigid piece. The script prints how far the drawn tip is from the
measured tip lift.

Output: results/figures/stage3_faces_meet.gif
"""

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import numpy as np
from PIL import Image

from analytical import Geometry, Material, effective_angles
from make_stage3_models import CHAMFER_DEPTH, TERM_DEPTH, terminal_geometry

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
OUT = RESULTS / "figures" / "stage3_faces_meet.gif"

MODEL = json.loads((RESULTS / "stage3_model.json").read_text())
X1, X2, X3, X4 = MODEL["tooth_stations"]
Y = MODEL["protrusion"]
L = MODEL["L"]
T = Geometry().t
ZA = T / 2                                   # beam axis, mid-thickness
UNDER = {int(k): v for k, v in MODEL["under_nodes"].items()}
XC = min(UNDER.values(), key=lambda x: abs(x - 0.5 * (X2 + X3)))   # crest node
RELEASE = Y - 1e-4                           # same test as extract_stage3.py
START = 1e-4                                 # lift that counts as contact

# case, direction, drawn angle, sign of the rotation, run label, face, colour
RUNS = (
    ("extract_mu020", +1, 45.0, -1, "Pull-out", "retention face", "#eb6834"),
    ("insert_mu020", -1, 30.0, +1, "Push-in", "lead-in face", "#2a78d6"),
)

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#ffffff"
LANCE_FILL, TERMINAL_FILL = "#e4e2db", "#c9c7bf"

SIZE = (9.6, 5.4)                            # inches at 100 dpi, so 960 x 540 px
VIEW_X, VIEW_Z = (3.0, 10.6), (-0.95, 1.65)   # mm; the root and the terminal's
                                             # lower part are outside the view
N_MOVE = 72                                  # frames from contact to release
MS_FIRST, MS_MOVE, MS_LAST = 1200, 60, 3600


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------

def load_runs():
    rows = {}
    with open(RESULTS / "stage3_curve_history.csv") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            rows.setdefault(r["case"], []).append(r)

    runs = []
    for case, direction, drawn, sign, label, face, colour in RUNS:
        rs = rows[case]
        lift = [abs(float(r["crest_lift_mm"])) for r in rs]
        rel = next(i for i, v in enumerate(lift) if v >= RELEASE)
        first = next(i for i, v in enumerate(lift) if v >= START)
        runs.append({
            "case": case, "direction": direction, "drawn": drawn, "sign": sign,
            "label": label, "face": face, "colour": colour,
            "start": max(0, first - 1), "release": rel,
            "lift": lift[:rel + 1],
            "rot": [abs(float(r["rotation_deg"])) for r in rs][:rel + 1],
            "travel": [float(r["travel_mm"]) for r in rs][:rel + 1],
            "tip": [abs(float(r["tip_lift_mm"])) for r in rs][:rel + 1],
        })
    return runs


def at_lift(run, h):
    """
    The run where its crest first reaches h, interpolated between the samples
    either side. Returns the state and the index of the later sample.
    """
    if h <= 0.0:
        i = run["start"]
        return {k: run[k][i] for k in ("lift", "rot", "travel", "tip")}, i + 1
    lift = run["lift"]
    k = next(i for i, v in enumerate(lift) if v >= h)
    f = (h - lift[k - 1]) / (lift[k] - lift[k - 1])
    state = {key: run[key][k - 1] + f * (run[key][k] - run[key][k - 1])
             for key in ("rot", "travel", "tip")}
    state["lift"] = h
    return state, k


def face_angle(run, rot):
    return run["drawn"] + run["sign"] * rot


def results_row(case):
    with open(RESULTS / "stage3_results.csv") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if r["case"] == case:
                return r
    return None


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

def line(a, b, n):
    return [(a[0] + (b[0] - a[0]) * i / n, a[1] + (b[1] - a[1]) * i / n)
            for i in range(n)]


def lance_outline():
    """Undeformed outline: top edge root to tip, then back along the underside."""
    pts = line((0.0, T), (L, T), 80)
    pts += line((L, T), (X4, 0.0), 4)
    pts += line((X4, 0.0), (X3, -Y), 10)       # lead-in
    pts += line((X3, -Y), (X2, -Y), 4)         # crest
    pts += line((X2, -Y), (X1, 0.0), 6)        # retention face
    pts += line((X1, 0.0), (0.0, 0.0), 60)
    pts.append((0.0, 0.0))
    return pts


def deform(points, lift, rot_deg):
    """
    Move points of the lance to the measured crest lift and rotation.

    Root to tooth: a cubic with zero lift and slope at the root, matching the
    lift and slope at x1. Tooth: a rigid turn about the axis at x1, with the
    lift at x1 chosen so the crest node lands on the measured crest lift.
    """
    th = math.radians(rot_deg)
    s, c = math.sin(th), math.cos(th)
    w1 = lift - (XC - X1) * s - (1.0 - c) * (ZA + Y)
    slope = math.tan(th)
    a = (3.0 * w1 - slope * X1) / X1 ** 2
    b = (slope * X1 - 2.0 * w1) / X1 ** 3
    out = []
    for x, z in points:
        if x <= X1:
            w = a * x ** 2 + b * x ** 3
            phi = math.atan(2.0 * a * x + 3.0 * b * x ** 2)
            out.append((x - (z - ZA) * math.sin(phi),
                        ZA + w + (z - ZA) * math.cos(phi)))
        else:
            dx, dz = x - X1, z - ZA
            out.append((X1 + dx * c - dz * s, ZA + w1 + dx * s + dz * c))
    return out


def terminal_outline(direction, travel):
    """The terminal's outline after it has moved, and the middle of its flat top."""
    g = terminal_geometry(direction)
    wall_a, wall_b = g["wall"]
    shift = direction * travel
    pts = [(g["lead"], -TERM_DEPTH), (g["lead"], -CHAMFER_DEPTH),
           (g["flat_a"], 0.0), (wall_a, 0.0), (wall_b, -Y),
           (g["end"], -Y), (g["end"], -TERM_DEPTH)]
    return [(x + shift, z) for x, z in pts], 0.5 * (g["flat_a"] + wall_a) + shift


def contact_face(run):
    if run["sign"] < 0:
        return [(X1, 0.0), (X2, -Y)]
    return [(X3, -Y), (X4, 0.0)]


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

def build_figure(runs):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
                         "xtick.color": MUTED, "ytick.color": MUTED,
                         "xtick.labelcolor": INK2, "ytick.labelcolor": INK2})
    fig = plt.figure(figsize=SIZE, dpi=100, facecolor=SURFACE)

    title = fig.text(0.030, 0.925, "", fontsize=15, color=INK, weight="bold",
                     va="baseline")
    fig.text(0.030, 0.870, "As the lance lifts it rotates. The same rotation makes "
             "the 45° face shallower and the 30° face steeper.",
             fontsize=10, color=INK2, va="baseline")

    art = {"title": title, "rows": []}
    outline0 = lance_outline()
    width = 0.500
    height = width * SIZE[0] / SIZE[1] * (VIEW_Z[1] - VIEW_Z[0]) / (VIEW_X[1] - VIEW_X[0])
    for run, bottom in zip(runs, (0.480, 0.105)):
        ax = fig.add_axes((0.030, bottom, width, height))
        ax.set_xlim(*VIEW_X)
        ax.set_ylim(*VIEW_Z)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")

        ghost = Polygon(outline0, closed=True, fill=False, edgecolor=MUTED,
                        lw=0.8, ls=(0, (3, 2)), zorder=2)
        ax.add_patch(ghost)
        term = Polygon([(0, 0)] * 3, closed=True, facecolor=TERMINAL_FILL,
                       edgecolor=INK2, lw=1.0, zorder=3)
        ax.add_patch(term)
        lance = Polygon(outline0, closed=True, facecolor=LANCE_FILL,
                        edgecolor=INK2, lw=1.0, zorder=4)
        ax.add_patch(lance)
        face, = ax.plot([], [], color=run["colour"], lw=3.4,
                        solid_capstyle="round", zorder=6)
        lance_txt = ax.text(0, 0, "lance", fontsize=8.5, color=INK2,
                            ha="center", va="center", zorder=7)
        arrow = "terminal →" if run["direction"] > 0 else "← terminal"
        term_txt = ax.text(0, -0.5, arrow, fontsize=8.5, color=INK,
                           ha="center", va="center", zorder=7, clip_on=True)

        # Row header: colour key on the left, text in ink.
        ax.add_line(Line2D([0.0, 0.045], [1.09, 1.09], transform=ax.transAxes,
                           color=run["colour"], lw=3.4, solid_capstyle="round",
                           clip_on=False))
        head = ax.text(0.062, 1.09, "", transform=ax.transAxes, fontsize=10.5,
                       color=INK, va="center")

        art["rows"].append({"ax": ax, "lance": lance, "term": term, "face": face,
                            "lance_txt": lance_txt, "term_txt": term_txt,
                            "head": head})

    ax = fig.add_axes((0.640, 0.170, 0.335, 0.575))
    ax.set_facecolor(SURFACE)
    ax.set_xlim(-0.015, 0.625)
    ax.set_ylim(28.0, 47.0)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    ax.set_xticklabels(["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6"])
    ax.set_yticks([30, 35, 40, 45])
    ax.set_yticklabels(["30°", "35°", "40°", "45°"])
    ax.tick_params(length=3, width=0.8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_xlabel("Lift at the crest (mm)", color=INK2, fontsize=9.5)
    ax.set_ylabel("Face angle from the insertion axis", color=INK2, fontsize=9.5)

    cursor = ax.axvline(0.0, color=AXIS, lw=1.0, zorder=1)
    ax.text(0.012, 45.55, "45° as drawn", fontsize=8.5, color=INK2, va="bottom")
    ax.text(0.012, 29.35, "30° as drawn", fontsize=8.5, color=INK2, va="top")

    curves = []
    for run in runs:
        ln, = ax.plot([], [], color=run["colour"], lw=2.0, solid_capstyle="round",
                      solid_joinstyle="round", zorder=4)
        dot, = ax.plot([], [], ls="none", marker="o", ms=9, mfc=run["colour"],
                       mec=SURFACE, mew=2.0, zorder=5)
        curves.append((ln, dot))

    handles = [Line2D([], [], color=r["colour"], lw=2.0,
                      label=f"{r['face'].capitalize()}, {r['label'].lower()} run")
               for r in runs]
    leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.02),
                    fontsize=8.5, frameon=False, borderaxespad=0.0,
                    handlelength=1.8, labelspacing=0.35)
    for t in leg.get_texts():
        t.set_color(INK2)

    final_head = ax.text(0.608, 32.75, "", fontsize=9, color=INK, ha="right",
                         va="center", weight="bold")
    final = ax.text(0.608, 31.1, "", fontsize=9, color=INK, ha="right",
                    va="center", linespacing=1.45)

    fig.text(0.030, 0.030,
             "Stage 3, µ = 0.20. Pull-out = extract_mu020, push-in = insert_mu020, "
             "from results/stage3_curve_history.csv.\n"
             "Face angle = drawn angle ∓ measured rotation. Lance drawn to scale "
             "from the measured lift and rotation, not the FE mesh. Runs matched by lift.",
             fontsize=7.6, color=MUTED, va="bottom", linespacing=1.4)

    art.update({"chart": ax, "cursor": cursor, "curves": curves,
                "final_head": final_head, "final": final})
    return fig, art


def draw_frame(art, runs, h, final_values=None):
    outline0 = lance_outline()
    worst_tip = 0.0
    for run, row, (ln, dot) in zip(runs, art["rows"], art["curves"]):
        state, k = at_lift(run, h)
        rot = state["rot"]
        pts = deform(outline0, state["lift"], rot)
        row["lance"].set_xy(pts)
        tip_drawn = deform([(X4, 0.0)], state["lift"], rot)[0][1]
        worst_tip = max(worst_tip, abs(tip_drawn - state["tip"]))

        term_pts, flat_mid = terminal_outline(run["direction"], state["travel"])
        row["term"].set_xy(term_pts)
        row["term_txt"].set_position((flat_mid, -0.5))

        fx, fz = zip(*deform(contact_face(run), state["lift"], rot))
        row["face"].set_data(fx, fz)
        lx, lz = deform([(VIEW_X[0] + 0.9, ZA)], state["lift"], rot)[0]
        row["lance_txt"].set_position((lx, lz))

        angle = face_angle(run, rot)
        row["head"].set_text(f"{run['label']}:  {run['face']}  {angle:.2f}°")

        xs = run["lift"][:k] + [state["lift"]]
        ys = [face_angle(run, r) for r in run["rot"][:k]] + [angle]
        if h <= 0.0:
            xs, ys = [state["lift"]], [angle]
        ln.set_data(xs, ys)
        dot.set_data([state["lift"]], [angle])

    art["cursor"].set_xdata([h, h])
    art["final_head"].set_text("" if final_values is None else "At release")
    art["final"].set_text("" if final_values is None else final_values)
    return worst_tip


def render(fig):
    fig.canvas.draw()
    return Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())


def to_palette(frames, colours=128):
    """
    One palette for every frame, so nothing flickers from frame to frame.

    The flat colours go in exactly, so white stays white and the two series
    colours are the ones that passed the palette check. The rest of the palette
    is fitted to the edge shades of three frames. Pixels are then mapped to the
    nearest entry by hand, because Pillow's own mapping is approximate and put
    the white background on a light grey.
    """
    keys = [tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))
            for h in (SURFACE, INK, INK2, MUTED, GRID, AXIS, LANCE_FILL,
                      TERMINAL_FILL, *(r[-1] for r in RUNS))]

    def packed(a):
        a = a.astype(np.int64)
        return (a[:, 0] << 16) | (a[:, 1] << 8) | a[:, 2]

    sample = np.concatenate([np.asarray(f).reshape(-1, 3)
                             for f in (frames[0], frames[len(frames) // 2], frames[-1])])
    rest = sample[~np.isin(packed(sample), packed(np.array(keys)))]
    rest = rest[:len(rest) // 1000 * 1000].reshape(-1, 1000, 3)
    n_fit = colours - len(keys)
    fitted = Image.fromarray(rest).quantize(colors=n_fit,
                                            method=Image.Quantize.MEDIANCUT)
    extra = fitted.getpalette()[:3 * n_fit]
    palette = np.array(keys + [tuple(extra[3 * i:3 * i + 3]) for i in range(n_fit)])

    out = []
    for f in frames:
        p = packed(np.asarray(f).reshape(-1, 3))
        uniq, inv = np.unique(p, return_inverse=True)
        u = np.stack([(uniq >> 16) & 255, (uniq >> 8) & 255, uniq & 255], axis=1)
        nearest = ((u[:, None, :] - palette[None, :, :]) ** 2).sum(axis=2).argmin(axis=1)
        im = Image.frombytes("P", f.size, nearest.astype(np.uint8)[inv].tobytes())
        im.putpalette(palette.astype(np.uint8).flatten().tolist())
        out.append(im)
    return out


def main():
    runs = load_runs()
    lead_in, retention = effective_angles(Geometry(), Material())
    predicted = (retention, lead_in)

    ends = [face_angle(r, at_lift(r, RELEASE)[0]["rot"]) for r in runs]
    for run, end in zip(runs, ends):
        row = results_row(run["case"])
        quoted = float(row["angle_from_shape_deg"]) if row else float("nan")
        flag = "" if abs(end - quoted) <= 0.01 else "   <- differs from the results file"
        print(f"{run['case']:<15} face at release {end:.2f} deg, "
              f"stage3_results.csv {quoted:.2f} deg{flag}")

    fig, art = build_figure(runs)
    art["title"].set_text(f"The two faces meet at {0.5 * (ends[0] + ends[1]):.1f}°")
    final_text = (f"measured  {ends[0]:.2f}° and {ends[1]:.2f}°\n"
                  f"predicted  {predicted[0]:.2f}° and {predicted[1]:.2f}°")

    lifts = [RELEASE * i / N_MOVE for i in range(N_MOVE + 1)]
    frames, durations, worst = [], [], 0.0
    for i, h in enumerate(lifts):
        last = i == len(lifts) - 1
        worst = max(worst, draw_frame(art, runs, h, final_text if last else None))
        frames.append(render(fig))
        durations.append(MS_FIRST if i == 0 else MS_LAST if last else MS_MOVE)
    plt.close(fig)

    frames = to_palette(frames)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(OUT, save_all=True, append_images=frames[1:], duration=durations,
                   loop=0, optimize=False, disposal=1)
    print(f"drawn tip against measured tip lift: worst {worst:.4f} mm")
    print(f"wrote {OUT.relative_to(REPO)}  ({len(frames)} frames, "
          f"{sum(durations) / 1000:.1f} s, {OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
