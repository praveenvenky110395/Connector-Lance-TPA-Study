#!/usr/bin/env python3
"""
Writes the Stage 3 LS-DYNA decks: the locking cycle.

Stage 3 stops being a beam exercise and becomes the part. The lance carries an
integral locking tooth - 45 deg retention face on the root side, flat crest,
30 deg lead-in towards the tip, protruding by y. The protrusion IS the lift
needed to release, which is the design relation the whole feature is built on. A
rigid terminal runs underneath with a matching shoulder, and a rigid TPA sits
above with a clearance.

Nothing new is assumed. The prediction comes from the two stages already done:

    how hard the lance is to lift          Stage 1, now on a variable section
    what the ramp does to the force        Stage 2, W/P = (tan b + mu)/(1 - mu tan b)

and three things change on the way from the idealisation to the part, all named
in analytical.py before any of this was run:

  1. the tooth stiffens the lance, by about 1 %
  2. the load sits on a tooth face, not at the tip, so the lever arm is shorter
     and the release condition belongs at the crest
  3. the inclined faces are on the lance, and the lance rotates

Item 3 is the finding. Stage 2 deliberately put the wedge face on the rigid body
so the contact normal could not rotate. On a real lance it is on the lance and it
does rotate - about 7.3 deg at release here - which takes the 30 deg lead-in to
roughly 37.8 and the 45 deg retention to roughly 37.7. The faces meet in the
middle and the asymmetry the two angles exist to produce largely disappears.

The nine runs

    lance_only        no terminal, no TPA. Prescribed tip lift. Must return the
                      Stage 1 FE force plus the ~1 % the tooth adds.
    extract_mu000     terminal pulled out, frictionless
    extract_mu020     same at mu = 0.20
    extract_mu030     same at mu = 0.30
    insert_mu020      terminal pushed in, mu = 0.20, full stroke through the snap
    extract_tpa_g010  extraction with the TPA at 0.10 mm clearance - blocked
    extract_tpa_g085  TPA at 0.85 mm, above the 0.757 mm blocking limit. Not
                      blocked, and it doubles as the control: a contact that is
                      defined but never reached has to cost exactly nothing.
    extract_mu020_fine  the baseline with the element over the tooth halved
    extract_mu020_slow  the baseline with the ramp doubled

Units: mm, tonne, s, N, MPa. Geometry is frozen in docs/numerical_spec.md.
"""

import json
import math
from pathlib import Path

from make_stage1_models import (
    L, B, T, Y_TIP, E, NU, RHO, NY,
    CONTROL_MASS, CONTROL_INERTIA, MASS_EID, INERTIA_EID,
)
from analytical import (Geometry, Material, tooth_stations, tooth_depth,
                        tpa_blocking_clearance, tpa_block_travel,
                        beam_flexibility, tip_force_timoshenko)

REPO = Path(__file__).resolve().parents[1]
G, M = Geometry(), Material()

N_THRU = 4                  # Stage 1's comparison mesh, through the thickness
TARGET_H = 0.20             # nominal element length along x
TOOTH_H = 0.10              # and over the tooth, where the contact lives. A
                            # master corner that is not detected penetrates by
                            # roughly the element length before a node notices,
                            # so halving it halves the worst case.

# Rigid parts. The modulus sizes the contact penalty stiffness and nothing else,
# so it is matched to the lance for the reason Stage 2 established: a steel value
# there would let the penalty spring rather than the structure set the time step.
RIGID_RHO = 7.85e-9
RIGID_E = E
RIGID_NU = 0.30

GAP_X = 0.010               # standoff between the two 45 deg faces at t = 0
TERM_FLAT = 1.60            # terminal flat top ahead of its shoulder
TERM_FLOOR = 1.60           # recess floor behind it
TERM_DEPTH = 1.60           # how far the terminal reaches below the lance
TERM_NZ = 3
OVERHANG = 0.25             # each side in y, so the rigid parts are wider

TPA_DEPTH = 0.60
TPA_NZ = 2
TPA_BACK = 0.30             # how far the TPA reaches behind the retention face
TPA_FRONT = 0.20            # and past the lance tip

EXTRACT_STROKE = 1.20       # The faces part at 0.537 mm of travel, not 0.600:
                            # the crest has to lift 0.600 and the contact sits
                            # nearer the root, where the lance lifts 0.894 as
                            # much. 1.20 leaves 0.66 mm after release so the
                            # tooth lands on the terminal and the force drops to
                            # friction only - the tail is what makes the peak
                            # readable as a peak.

# How far past the stop a blocked run is driven.
#
# This is the whole reason the first TPA runs failed. The terminal is rigid and
# kinematically driven: once the lance is resting on the TPA there is nothing in
# the model that can stop the terminal, so it simply ploughs on, drags the tooth
# along x and stretches the lance. The first attempt drove 0.90 mm against a
# 0.079 mm stop - eleven times past it.
#
# A blocked run is therefore driven to the stop plus this much and no further.
# The overrun is not padding: it is where the force rises steeply and the load
# path moves off the lance's bending and into the TPA, which is the thing the
# stage set out to measure. Beyond it the linear-elastic lance has fractured and
# the model is describing a part that no longer exists.
BLOCK_OVERRUN = 0.05
INSERT_LEAD = 0.10          # where the terminal corner starts past the lance tip

# Lead-in chamfer on the terminal's leading edge. A real terminal has one, and
# numerically it turns a single 90 deg corner into two 135 deg ones, which is a
# far gentler thing to drag along a deformable face. Kept steep and short on
# purpose: the tooth's own 30 deg lead-in still has to be the shallower of the
# two surfaces, because that is the angle the insertion force is predicted from.
CHAMFER_DEPTH = 0.15
CHAMFER_ANGLE = 45.0

T_RAMP = 5.0e-3
T_END = 5.5e-3
T_RAMP_INSERT = 30.0e-3     # longer, because the insertion stroke is four times
T_END_INSERT = 31.0e-3      # as long and the speed is what the energy check sees.
                            # Raised from 20 ms after the first attempt: peak
                            # speed there was 352 mm/s against Stage 2's 212,
                            # and a contact that is already marginal is worse at
                            # speed.

PID_LANCE, PID_TERMINAL, PID_TPA, PID_TIP_BODY = 1, 2, 3, 4
TERM_NODE0, TERM_ELEM0 = 3000, 3000
TPA_NODE0, TPA_ELEM0 = 6000, 6000

# (name, mu, tpa clearance or None, direction, options)
#   direction  +1 extract, -1 insert, 0 no terminal
#   options    tooth_h overrides the element size over the tooth
#              t_ramp  overrides the ramp time
#
# The last two runs are the sensitivity pair, and they exist because a contact
# result that has not been shown to be independent of the mesh and of the rate
# is not a result. Both are extract_mu020 with exactly one thing changed:
#
#   _fine  halves the element length over the tooth. A penalty contact
#          transmits force by penetrating, and the penetration scales with the
#          element, so this is the direct test of whether the force is set by
#          the mechanics or by the discretisation.
#   _slow  doubles the ramp time. The energy check says the run is quasi-static;
#          this is the test that does not take the energy check's word for it.
#
# Neither is expected to move the answer. Running them is what lets that be
# said rather than assumed.
CASES = [
    ("lance_only",       0.00, None,  0, {}),
    ("extract_mu000",    0.00, None, +1, {}),
    ("extract_mu020",    0.20, None, +1, {}),
    ("extract_mu030",    0.30, None, +1, {}),
    ("insert_mu020",     0.20, None, -1, {}),
    ("extract_tpa_g010", 0.20, 0.10, +1, {}),
    # 0.85, not 0.65. The TPA reaches the lance tip, and the tip lifts 1.26
    # times the crest, so it blocks below 0.757 mm of clearance rather than
    # below the 0.60 mm protrusion. The 0.65 mm run that was meant to show a
    # fitted and useless TPA in fact blocked - see tpa_blocking_clearance.
    ("extract_tpa_g085", 0.20, 0.85, +1, {}),
    ("extract_mu020_fine", 0.20, None, +1, {"tooth_h": TOOTH_H / 2.0}),
    ("extract_mu020_slow", 0.20, None, +1, {"t_ramp": 2.0 * T_RAMP}),
]


# ---------------------------------------------------------------------------
# Lance with the tooth
# ---------------------------------------------------------------------------

def x_stations(tooth_h=None):
    """
    Element boundaries along the lance, with the four tooth breakpoints landing
    exactly on one. The faces are then planes rather than staircases, which is
    the whole reason the mesh is graded instead of uniform.
    """
    th = TOOTH_H if tooth_h is None else tooth_h
    x1, x2, x3, x4 = tooth_stations(G)
    out = []
    for a, b, h in ((0.0, x1, TARGET_H), (x1, x2, th),
                    (x2, x3, th), (x3, x4, th)):
        n = max(1, round((b - a) / h))
        out += [a + (b - a) * i / n for i in range(n)]
    out.append(x4)
    return out


def build_lance(tooth_h=None):
    """
    Structured hex mesh whose underside follows the tooth profile. The top stays
    flat at T, so the section thickens under the tooth - which is the part, and
    is also why the lance comes out about 1 % stiffer than the Stage 1 prism.
    """
    xs = x_stations(tooth_h)
    nx, nz = len(xs) - 1, N_THRU

    def nid(i, j, k):
        return i * (NY + 1) * (nz + 1) + j * (nz + 1) + k + 1

    nodes = []
    for i, x in enumerate(xs):
        z_bot = -tooth_depth(G, x)
        for j in range(NY + 1):
            for k in range(nz + 1):
                nodes.append((nid(i, j, k), x, B * j / NY,
                              z_bot + (T - z_bot) * k / nz))

    elements = []
    for i in range(nx):
        for j in range(NY):
            for k in range(nz):
                elements.append((len(elements) + 1, [
                    nid(i, j, k), nid(i + 1, j, k),
                    nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                ]))

    control = len(nodes) + 1
    nodes.append((control, L, B / 2.0, T / 2.0))

    root = sorted(nid(0, j, k) for j in range(NY + 1) for k in range(nz + 1))
    tip = sorted(nid(nx, j, k) for j in range(NY + 1) for k in range(nz + 1))
    spine = [nid(i, NY // 2, nz) for i in range(nx + 1)]
    under = [nid(i, NY // 2, 0) for i in range(nx + 1)]

    return {"xs": xs, "nx": nx, "nz": nz, "nodes": nodes, "elements": elements,
            "root": root, "tip": tip, "spine": spine, "under": under,
            "control": control}


# ---------------------------------------------------------------------------
# Rigid parts
# ---------------------------------------------------------------------------

def profile_block(xs, top_of, z_bottom, ny, nz, node0, elem0):
    """Structured block under a given top profile."""
    def nid(i, j, k):
        return node0 + i * (ny + 1) * (nz + 1) + j * (nz + 1) + k

    y0, y1 = -OVERHANG, B + OVERHANG
    nodes = []
    for i, x in enumerate(xs):
        z_top = top_of(x)
        for j in range(ny + 1):
            for k in range(nz + 1):
                nodes.append((nid(i, j, k), x, y0 + (y1 - y0) * j / ny,
                              z_bottom + (z_top - z_bottom) * k / nz))

    elements = []
    for i in range(len(xs) - 1):
        for j in range(ny):
            for k in range(nz):
                elements.append((elem0 + len(elements), [
                    nid(i, j, k), nid(i + 1, j, k),
                    nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                ]))
    return nodes, elements


def extraction_stroke(gap):
    """
    How far to pull the terminal. Full stroke unless the TPA stops the lance
    first, in which case the stop plus BLOCK_OVERRUN and not a micron more.
    """
    if gap is None or not _blocks(gap):
        return EXTRACT_STROKE
    return tpa_block_travel(G, M, gap, tpa_reach()) + BLOCK_OVERRUN


def tpa_reach():
    """The furthest station along the lance that the TPA covers."""
    _, _, _, x4 = tooth_stations(G)
    return min(L, x4 + TPA_FRONT)


def _blocks(gap):
    return gap < tpa_blocking_clearance(G, M, tpa_reach())


def terminal_geometry(direction, gap=None):
    """
    The terminal, in the position it starts the run in.

    Its top is flat at the lance underside, then a 45 deg shoulder drops to a
    recess the tooth sits in. The shoulder is the same plane as the tooth's
    retention face, offset by GAP_X so the run does not open on an initial
    penetration.

    Extraction starts seated and pulls +x. Insertion starts with the leading
    corner just past the lance tip and pushes -x until the tooth drops into the
    recess, so the whole ride up the lead-in, the slide along the crest and the
    snap are all in one stroke.
    """
    x1, x2, _, x4 = tooth_stations(G)
    seated_lead = x1 - GAP_X - TERM_FLAT          # leading corner when seated
    chamfer_run = CHAMFER_DEPTH / math.tan(math.radians(CHAMFER_ANGLE))
    seated_lead -= chamfer_run                    # the chamfer sits ahead of the flat
    if direction < 0:
        lead = x4 + INSERT_LEAD
        stroke = -(lead - seated_lead)
    else:
        lead = seated_lead
        stroke = extraction_stroke(gap)

    chamfer = CHAMFER_DEPTH / math.tan(math.radians(CHAMFER_ANGLE))
    flat_a = lead + chamfer                       # chamfer meets the flat top
    wall_a = flat_a + TERM_FLAT                   # top of the 45 deg shoulder
    wall_b = wall_a + G.y / math.tan(math.radians(G.beta_deg))
    end = wall_b + TERM_FLOOR

    def top_of(x):
        if x <= flat_a:
            return -CHAMFER_DEPTH + (x - lead) * math.tan(math.radians(CHAMFER_ANGLE))
        if x <= wall_a:
            return 0.0
        if x >= wall_b:
            return -G.y
        return -(x - wall_a) * math.tan(math.radians(G.beta_deg))

    xs = []
    for a, b in ((lead, flat_a), (flat_a, wall_a), (wall_a, wall_b), (wall_b, end)):
        n = max(1, round((b - a) / TARGET_H))
        xs += [a + (b - a) * i / n for i in range(n)]
    xs.append(end)

    nodes, elements = profile_block(xs, top_of, -TERM_DEPTH, NY, TERM_NZ,
                                    TERM_NODE0, TERM_ELEM0)
    return {"nodes": nodes, "elements": elements, "stroke": stroke,
            "wall": (wall_a, wall_b), "lead": lead, "end": end,
            "flat_a": flat_a, "top_of": top_of}


def tpa_geometry(clearance):
    """Rigid block above the lance, flat underside at the stated clearance."""
    x1, _, _, x4 = tooth_stations(G)
    x_a, x_b = x1 - TPA_BACK, x4 + TPA_FRONT
    n = max(1, round((x_b - x_a) / TARGET_H))
    xs = [x_a + (x_b - x_a) * i / n for i in range(n + 1)]
    nodes, elements = profile_block(
        xs, lambda x: T + clearance + TPA_DEPTH, T + clearance,
        NY, TPA_NZ, TPA_NODE0, TPA_ELEM0)
    return {"nodes": nodes, "elements": elements, "underside": T + clearance}


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_volume(nodes, elements, expected, what):
    coords = {n[0]: n[1:] for n in nodes}

    def tet(a, b, c, d):
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        w = [d[i] - a[i] for i in range(3)]
        return (u[0] * (v[1] * w[2] - v[2] * w[1])
                - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0])) / 6.0

    tets = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6), (3, 4, 6, 7), (1, 4, 5, 6)]
    total = 0.0
    for _, conn in elements:
        p = [coords[n] for n in conn]
        v = sum(tet(p[a], p[b], p[c], p[d]) for a, b, c, d in tets)
        if v <= 0.0:
            raise ValueError(f"negative Jacobian in the {what} mesh")
        total += v
    if expected is not None and abs(total - expected) > 1e-6 * expected:
        raise ValueError(f"{what} volume {total:.6f}, expected {expected:.6f}")
    return total


def check_lance(lance):
    """The tooth has to come out of the mesh at the angles it was drawn at."""
    x1, x2, x3, x4 = tooth_stations(G)
    coords = {n[0]: n[1:] for n in lance["nodes"][:-1]}
    for station in (x1, x2, x3, x4):
        if not any(abs(x - station) < 1e-9 for x in lance["xs"]):
            raise ValueError(f"tooth station {station:.4f} is not a mesh station")

    under = [coords[n] for n in lance["under"]]
    for a, b, ang in ((x1, x2, G.beta_deg), (x3, x4, G.alpha_deg)):
        pts = [(p[0], p[2]) for p in under if a - 1e-9 <= p[0] <= b + 1e-9]
        if len(pts) < 2:
            raise ValueError(f"face {a:.3f}-{b:.3f} has too few mesh points")
        for (xa, za), (xb, zb) in zip(pts, pts[1:]):
            slope = abs((zb - za) / (xb - xa))
            if abs(slope - math.tan(math.radians(ang))) > 1e-9:
                raise ValueError(f"face {a:.3f}-{b:.3f} is not at {ang} deg")

    crest = [p[2] for p in under if x2 + 1e-9 < p[0] < x3 - 1e-9]
    if crest and max(abs(z + G.y) for z in crest) > 1e-9:
        raise ValueError("crest is not at the full protrusion")

    # Volume: prism plus the tooth's own cross-section swept across the width.
    tooth_area = (0.5 * (x2 - x1) * G.y + (x3 - x2) * G.y
                  + 0.5 * (x4 - x3) * G.y)
    check_volume(lance["nodes"][:-1], lance["elements"],
                 L * B * T + tooth_area * B, "lance")


def check_clearances(lance, terminal, tpa, direction):
    """Nothing may start inside anything else."""
    coords = {n[0]: n[1:] for n in lance["nodes"][:-1]}
    lance_pts = list(coords.values())

    if terminal is not None:
        wall_a, wall_b = terminal["wall"]
        if direction > 0:
            # Seated: the shoulder plane must sit GAP_X ahead of the tooth face,
            # measured along x, and must not have crossed it.
            x1, _, _, _ = tooth_stations(G)
            if abs((x1 - GAP_X) - wall_a) > 1e-9:
                raise ValueError(f"shoulder at {wall_a:.4f}, expected "
                                 f"{x1 - GAP_X:.4f}")
        else:
            if terminal["lead"] <= L:
                raise ValueError("insertion starts with the corner inside the tooth")
        # No lance node may start below the terminal's top profile, evaluated at
        # that node's own x. An earlier version of this compared against the
        # highest terminal node within one element length, which is fine over a
        # flat top and wrong over a 45 deg wall - it read the flat ahead of the
        # wall as the local surface and failed the seated position, where every
        # node is correctly one GAP_X clear.
        top_of = terminal["top_of"]
        for x, _y, z in lance_pts:
            if not (terminal["lead"] - 1e-9 <= x <= terminal["end"] + 1e-9):
                continue
            if z < top_of(x) - 1e-9:
                raise ValueError(f"lance node at x={x:.4f}, z={z:.4f} starts "
                                 f"{top_of(x) - z:.4f} mm inside the terminal, "
                                 f"whose top there is {top_of(x):.4f}")

    if tpa is not None:
        gap = tpa["underside"] - T
        limit = tpa_blocking_clearance(G, M, tpa_reach())
        if abs(gap / limit - 1.0) < 0.10:
            raise ValueError(
                f"TPA clearance {gap:.3f} is within 10 % of the blocking limit "
                f"{limit:.3f} - the run would not clearly show either outcome")
        if tpa["underside"] <= T:
            raise ValueError("TPA underside is at or below the lance top face")
        if max(z for _x, _y, z in lance_pts) > tpa["underside"] - 1e-9:
            raise ValueError("lance starts inside the TPA")


def check_ids(text, expect_parts):
    """
    Read the finished deck back and check what its IDs point at.

    This exists because the previous Stage 3 decks were card-perfect and did
    nothing: the tip rigid body and a rigid part had both been given part 2, so
    the prescribed motion landed on a part that MAT_RIGID had locked in all six
    directions. LS-DYNA ran it to termination and wrote zeros. Counting cards
    could not see it; only checking what the IDs refer to can.
    """
    lines = [ln for ln in text.splitlines() if not ln.startswith("$")]

    def fields(i, w=10):
        return [lines[i][j:j + w].strip() for j in range(0, len(lines[i]), w)]

    # *MAT_RIGID con1 codes, as sets of the directions each one locks.
    con1_locks = {0: set(), 1: {1}, 2: {2}, 3: {3}, 4: {1, 2},
                  5: {2, 3}, 6: {3, 1}, 7: {1, 2, 3}}

    parts, sections, materials, sets = [], [], [], []
    bodies, prescribed, contacts, curves = [], [], [], []
    rigid_con1 = {}
    part_material = {}
    for i, ln in enumerate(lines):
        if ln.startswith("*PART"):
            row = fields(i + 2)
            parts.append(int(row[0]))
            part_material[int(row[0])] = int(row[2])
        elif ln.startswith("*SECTION_SOLID"):
            sections.append(int(fields(i + 1)[0]))
        elif ln.startswith("*MAT_ELASTIC"):
            materials.append(int(fields(i + 1)[0]))
        elif ln.startswith("*MAT_RIGID"):
            mid = int(fields(i + 1)[0])
            materials.append(mid)
            rigid_con1[mid] = int(fields(i + 2)[1])
        elif ln.startswith("*SET_NODE_LIST"):
            sets.append(int(fields(i + 1)[0]))
        elif ln.startswith("*CONSTRAINED_NODAL_RIGID_BODY"):
            bodies.append(int(fields(i + 1)[0]))
        elif ln.startswith("*BOUNDARY_PRESCRIBED_MOTION_RIGID"):
            prescribed.append([int(v) for v in fields(i + 1)[:4]])
        elif ln.startswith("*CONTACT_"):
            contacts.append([int(v) for v in fields(i + 1)[:4]])
        elif ln.startswith("*DEFINE_CURVE"):
            curves.append(int(fields(i + 1)[0]))

    if len(parts + bodies) != len(set(parts + bodies)):
        raise ValueError(f"part IDs collide: parts {parts}, bodies {bodies}")
    for name, ids in (("section", sections), ("material", materials),
                      ("node set", sets), ("curve", curves)):
        if len(ids) != len(set(ids)):
            raise ValueError(f"{name} IDs are not unique: {ids}")
    for pid, _dof, _vad, lcid in prescribed:
        if pid not in parts + bodies:
            raise ValueError(f"prescribed motion on unknown part {pid}")
        if pid in parts and pid == PID_LANCE:
            raise ValueError("prescribed motion on the deformable lance")
        if lcid not in curves:
            raise ValueError(f"prescribed motion uses curve {lcid}, not defined")
        # The one that cost a run: a rigid part driven along a direction its own
        # material card has already constrained. LS-DYNA does not complain - it
        # just does not move.
        con1 = rigid_con1.get(part_material.get(pid))
        if con1 is not None and _dof in con1_locks.get(con1, set()):
            raise ValueError(
                f"part {pid} is driven along dof {_dof} but *MAT_RIGID con1 = "
                f"{con1} constrains {sorted(con1_locks[con1])}")
    for ssid, msid, sstyp, mstyp in contacts:
        for sid, styp in ((ssid, sstyp), (msid, mstyp)):
            if styp == 3 and sid not in parts:
                raise ValueError(f"contact names part {sid}, which is not a *PART")
    if sorted(parts) != sorted(expect_parts):
        raise ValueError(f"expected parts {expect_parts}, deck has {parts}")


# ---------------------------------------------------------------------------
# Deck
# ---------------------------------------------------------------------------

def smoothstep_curve(t_ramp, t_hold, value, n=501):
    """Quintic smoothstep to `value`, then held past termination."""
    pts = []
    for i in range(n):
        s = i / (n - 1)
        pts.append((s * t_ramp, value * (10 * s**3 - 15 * s**4 + 6 * s**5)))
    pts.append((t_hold, value))
    return pts


def write_deck(lance, terminal, tpa, name, mu, direction, path, opts=None):
    opts = opts or {}
    insert = direction < 0
    t_ramp = T_RAMP_INSERT if insert else T_RAMP
    t_end = T_END_INSERT if insert else T_END
    if "t_ramp" in opts:
        t_ramp = opts["t_ramp"]
        t_end = t_ramp * (T_END / T_RAMP)
    # Two output rates, and the reason is the first Stage 3 run.
    #
    # Everything was written at 1e-5 s. The contact force rings far above
    # 50 kHz, so at that rate it was aliased: the apparent frequency came out
    # 44 kHz in one run and 16 kHz in the same model run slower, which a real
    # structural mode cannot do. Folded noise cannot be filtered out
    # afterwards, and every force number from that set had to be thrown away.
    #
    # Forces are now written at 1 MHz. Displacements stay at the old rate -
    # they were never the problem, because a displacement is the double
    # integral of the acceleration and the ring barely shows in it. That is
    # also why the rotation measurement survived the first set intact.
    state = 4.0e-5 if insert else 1.0e-5
    force = 5.0e-7
    # One rate for the forces on every run, insertion included. It was first
    # set slower there, on the reasoning that the insertion run is six times
    # longer and loaded six times slower. The resolution check in the extractor
    # rejected it: the ring is a property of the contact and the structure, not
    # of the loading rate, so halving the sample rate halves the resolution and
    # nothing else.
    #
    # 0.5 us, not 1 us. At 1 us the ring came back at 9 to 12 samples per cycle
    # on some runs and 7 on others, which is at the edge; at 0.5 us it is 12 to
    # 19 everywhere. That costs about 30 MB per extraction run and 90 MB for
    # the insertion one, which is a fair price for not having to wonder.

    parts = [PID_LANCE]
    lines = [
        f"$ Stage 3 locking cycle - {name}",
        f"$ mu {mu:.2f}, "
        + ("no terminal" if terminal is None else
           ("insertion" if insert else "extraction"))
        + (", no TPA" if tpa is None else
           f", TPA clearance {tpa['underside'] - T:.3f} mm"),
        "*KEYWORD",
        "*TITLE",
        f"Stage3 {name}",
        "*CONTROL_TERMINATION",
        f"{t_end:10.4E}",
        "*CONTROL_TIMESTEP",
        "$#  dtinit    tssfac      isdo    tslimt     dt2ms",
        "       0.0       0.9         0       0.0       0.0",
        "*CONTROL_ENERGY",
        "$#    hgen      rwen    slnten     rylen",
        "         2         2         2         2",
    ]
    for db in ("GLSTAT", "MATSUM", "NODOUT", "BNDOUT"):
        lines += [f"*DATABASE_{db}", f"{state:10.4E}"]
    for db in ("SPCFORC", "RCFORC"):
        lines += [f"*DATABASE_{db}", f"{force:10.4E}"]
    lines += ["*DATABASE_BINARY_D3PLOT", f"{t_ramp / 40.0:10.4E}"]

    lines += [
        "*PART", "Lance",
        "$#     pid     secid       mid     eosid      hgid",
        f"{PID_LANCE:10d}         1         1         0         0",
        "*SECTION_SOLID", "         1        -1",
        "*MAT_ELASTIC",
        "$#     mid        ro         e        pr",
        f"         1{RHO:10.3E}{E:10.1f}{NU:10.2f}",
    ]

    def rigid_part(pid, label, secid, mid, con1):
        """
        con1 is *MAT_RIGID's translational constraint code with cmo = 1.0, and
        it is a code, not a bitmask - which is how the terminal ended up locked.

            0 none   1 x   2 y   3 z   4 x+y   5 y+z   6 z+x   7 x+y+z

        The terminal is driven along x, so it needs 5, leaving x free. It was
        first written as 6, read as "everything except x". 6 constrains z and x,
        so the prescribed motion had nothing left to move: the run finished with
        the terminal exactly where it started, no contact, no force. Stage 2's
        plate is driven the same way and had 5 in it the whole time.
        """
        return [
            "*PART", label,
            "$#     pid     secid       mid     eosid      hgid",
            f"{pid:10d}{secid:10d}{mid:10d}         0         0",
            "*SECTION_SOLID", f"{secid:10d}         1",
            "*MAT_RIGID",
            "$#     mid        ro         e        pr",
            f"{mid:10d}{RIGID_RHO:10.3E}{RIGID_E:10.1f}{RIGID_NU:10.2f}",
            "$#     cmo      con1      con2",
            f"       1.0{con1:10d}         7",
            "$ third card, mandatory even when empty",
            "         0       0.0       0.0       0.0       0.0       0.0       0.0",
        ]

    if terminal is not None:
        # con1 = 6 leaves x free, which the prescribed motion then drives.
        # 5 = y and z constrained, x free for the prescribed motion to drive.
        lines += rigid_part(PID_TERMINAL, "Terminal", 2, 2, 5)
        parts.append(PID_TERMINAL)
    if tpa is not None:
        lines += rigid_part(PID_TPA, "TPA", 3, 3, 7)
        parts.append(PID_TPA)

    lines.append("*NODE")
    for group in (lance["nodes"],
                  terminal["nodes"] if terminal else [],
                  tpa["nodes"] if tpa else []):
        for n, x, y, z in group:
            lines.append(f"{n:8d}{x:16.8f}{y:16.8f}{z:16.8f}")

    lines.append("*ELEMENT_SOLID")
    for group, pid in ((lance["elements"], PID_LANCE),
                       (terminal["elements"] if terminal else [], PID_TERMINAL),
                       (tpa["elements"] if tpa else [], PID_TPA)):
        for eid, conn in group:
            lines.append(f"{eid:8d}{pid:8d}")
            lines.append("".join(f"{n:8d}" for n in conn))

    sets = [(1, lance["root"]), (2, lance["tip"]), (3, lance["spine"]),
            (4, [lance["control"]] + lance["tip"]), (5, lance["under"])]
    for sid, ids in sets:
        lines += ["*SET_NODE_LIST", f"{sid:10d}"]
        for i in range(0, len(ids), 8):
            lines.append("".join(f"{n:10d}" for n in ids[i:i + 8]))

    lines += [
        "*DATABASE_HISTORY_NODE_SET",
        "$ Set 3 is the top surface root to tip, set 5 the underside - the",
        "$ underside is where the tooth is, so that is the line that says how",
        "$ far the crest has lifted and whether it has cleared.",
        "         3         5",
        "*DATABASE_HISTORY_NODE",
        f"{lance['control']:10d}",
        "*BOUNDARY_SPC_SET",
        "$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz",
        "         1         0         1         1         1         1         1         1",
    ]

    if terminal is None:
        # The Stage 1 driver, lifting the tip. This run exists to show that the
        # tooth has not changed the beam into something else.
        lines += [
            "*CONSTRAINED_NODAL_RIGID_BODY",
            "$#     pid       cid      nsid     pnode      iprt    drflag    rrflag",
            f"{PID_TIP_BODY:10d}         0         4         0         0         0         0",
            "*ELEMENT_MASS",
            f"{MASS_EID:8d}{lance['control']:8d}{CONTROL_MASS:16.4E}",
            "*ELEMENT_INERTIA",
            f"{INERTIA_EID:8d}{lance['control']:8d}",
            f"{CONTROL_INERTIA:10.3E}{0.0:10.1f}{0.0:10.1f}"
            f"{CONTROL_INERTIA:10.3E}{0.0:10.1f}{CONTROL_INERTIA:10.3E}",
            "*BOUNDARY_PRESCRIBED_MOTION_RIGID",
            "$#     pid       dof       vad      lcid        sf",
            f"{PID_TIP_BODY:10d}         3         2         1       1.0",
        ]
        curve = smoothstep_curve(t_ramp, t_end + 1.0e-3, Y_TIP)
    else:
        lines += [
            "*BOUNDARY_PRESCRIBED_MOTION_RIGID",
            "$ The terminal is driven along its own axis. Nothing is prescribed",
            "$ on the lance at all - it is loaded only through the tooth.",
            "$#     pid       dof       vad      lcid        sf",
            f"{PID_TERMINAL:10d}         1         2         1       1.0",
            "*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE",
            "$ TWO-WAY, and the one-way version is why the first insertion run",
            "$ died of element distortion.",
            "$",
            "$ One-way checks slave NODES against master SEGMENTS and nothing",
            "$ else. The terminal's leading edge is a sharp rigid corner that",
            "$ travels the whole length of the tooth, and a master corner can sit",
            "$ inside a slave element face, between its nodes, entirely",
            "$ undetected. It gouged its way along the lead-in face from about",
            "$ 4 ms and the run terminated at 10.1 ms with the corner under the",
            "$ retention face. Two-way checks both sides, so the corner is a node",
            "$ against the lance's faces as well.",
            "$",
            "$ Stage 2 put the ramp on the rigid master so the normal could not",
            "$ rotate. Here it is on the lance, because that is where it is on the",
            "$ part - the rotation of the normal is the thing being measured.",
            "$#    ssid      msid     sstyp     mstyp",
            f"{PID_LANCE:10d}{PID_TERMINAL:10d}         3         3",
            "$#      fs        fd        dc        vc       vdc    penchk",
            f"{mu:10.2f}{mu:10.2f}       0.0       0.0       0.0         0",
            "$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf",
            "       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0",
        ]
        curve = smoothstep_curve(t_ramp, t_end + 1.0e-3, terminal["stroke"])

    if tpa is not None:
        lines += [
            "*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE",
            "$ Lance against the TPA, two-way for the same reason.",
            "$ Frictionless: the TPA is a stop, and",
            "$ friction there would put an axial force into a part whose only",
            "$ job is to be in the way.",
            "$#    ssid      msid     sstyp     mstyp",
            f"{PID_LANCE:10d}{PID_TPA:10d}         3         3",
            "$#      fs        fd        dc        vc       vdc    penchk",
            "      0.00      0.00       0.0       0.0       0.0         0",
            "$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf",
            "       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0",
        ]

    lines += ["*DEFINE_CURVE", "         1"]
    for t, d in curve:
        lines.append(f"{t:20.10E}{d:20.10E}")
    lines.append("*END")

    text = "\n".join(lines) + "\n"
    check_ids(text, parts)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return len(text.splitlines())


def stage1_reference():
    """
    Stage 1's measured force at the mesh Stage 3 uses, from its results file.

    Stage 3's lance has N_THRU elements through the thickness, so the row that
    applies is the Stage 1 mesh run with the same count. Returns None if Stage
    1 has not been evaluated, and the caller then falls back to beam theory
    and says so rather than inventing a baseline.
    """
    path = REPO / "results" / "stage1_results.csv"
    if not path.exists():
        return None
    import csv
    with open(path) as f:
        for row in csv.DictReader(l for l in f if not l.startswith("#")):
            if (row.get("n_thru") == str(N_THRU) and row.get("elform") == "-1"
                    and row.get("peak_force_N", "").strip()):
                return float(row["peak_force_N"])
    return None


def write_run_scripts(out, names):
    """
    Two loops over the case folders, one per shell. Each one changes into the
    folder before it starts the solver, because that is where LS-DYNA puts its
    output and the point of the folders is that the nine runs do not overwrite
    each other.
    """
    bat = [
        "@echo off",
        "rem Runs every Stage 3 deck in its own folder.",
        "rem   run_all.bat \"C:\\path\\to\\lsdyna.exe\" 4",
        "setlocal",
        "set SOLVER=%~1",
        "if \"%SOLVER%\"==\"\" set SOLVER=lsdyna",
        "set NCPU=%~2",
        "if \"%NCPU%\"==\"\" set NCPU=4",
        "",
    ]
    for n in names:
        bat += [
            f"echo === stage3_{n} ===",
            f"pushd \"%~dp0stage3_{n}\"",
            f"\"%SOLVER%\" i=stage3_{n}.k ncpu=%NCPU% memory=200m",
            "popd",
            "",
        ]
    bat += ["echo Done. Next:",
            "echo   python scripts\\extract_stage3.py \"%~dp0\"",
            "echo   python scripts\\postprocess_stage3.py"]
    (out / "run_all.bat").write_text("\r\n".join(bat) + "\r\n")

    sh = [
        "#!/bin/sh",
        "# Runs every Stage 3 deck in its own folder.",
        "#   ./run_all.sh /path/to/ls-dyna 4",
        "set -e",
        "SOLVER=\"${1:-ls-dyna}\"",
        "NCPU=\"${2:-4}\"",
        "HERE=\"$(cd \"$(dirname \"$0\")\" && pwd)\"",
        "",
    ]
    for n in names:
        sh += [
            f"echo \"=== stage3_{n} ===\"",
            f"(cd \"$HERE/stage3_{n}\" && \"$SOLVER\" i=stage3_{n}.k ncpu=\"$NCPU\" memory=200m)",
            "",
        ]
    sh += ["echo 'Done. Next:'",
           "echo '  python scripts/extract_stage3.py \"'\"$HERE\"'\"'",
           "echo '  python scripts/postprocess_stage3.py'"]
    p = out / "run_all.sh"
    p.write_text("\n".join(sh) + "\n")
    p.chmod(0o755)


def main():
    lance = build_lance()
    check_lance(lance)
    x1, x2, x3, x4 = tooth_stations(G)
    print(f"Lance: {len(lance['nodes'])} nodes, {len(lance['elements'])} elements, "
          f"{lance['nx']} along x")
    print(f"Tooth: retention {x1:.4f}-{x2:.4f} at {G.beta_deg:.0f} deg, "
          f"crest {x2:.4f}-{x3:.4f}, lead-in {x3:.4f}-{x4:.4f} at "
          f"{G.alpha_deg:.0f} deg\n")

    out = REPO / "ls-dyna" / "stage3"

    def node_map(lc):
        c = {n[0]: n[1:] for n in lc["nodes"]}
        return {"root_nodes": {str(n): c[n][2] for n in lc["root"]},
                "under_nodes": {str(n): c[n][0] for n in lc["under"]},
                "spine_nodes": {str(n): c[n][0] for n in lc["spine"]}}

    index = {"control_node": lance["control"], "L": L, "protrusion": G.y,
             "tooth_stations": list(tooth_stations(G)),
             # The extractor needs this to know where the terminal's nose stops
             # carrying the contact and the tooth's own lead-in takes over.
             "chamfer_depth": CHAMFER_DEPTH,
             "cases": {}}
    index.update(node_map(lance))

    print(f"{'case':>20} {'mu':>5} {'tpa gap':>8} {'stroke':>8} {'tooth h':>8} "
          f"{'ramp ms':>8} {'elements':>9}")
    for name, mu, gap, direction, opts in CASES:
        # A run that changes the mesh needs its own lance; the rest share one.
        this_lance = build_lance(opts["tooth_h"]) if "tooth_h" in opts else lance
        if this_lance is not lance:
            check_lance(this_lance)
            this_lance["under"] = this_lance["under"]

        terminal = terminal_geometry(direction, gap) if direction else None
        tpa = tpa_geometry(gap) if gap is not None else None
        check_clearances(this_lance, terminal, tpa, direction)
        if terminal:
            check_volume(terminal["nodes"], terminal["elements"], None, "terminal")
        if tpa:
            check_volume(tpa["nodes"], tpa["elements"], None, "TPA")

        # One folder per case. LS-DYNA writes d3plot, glstat, spcforc, nodout
        # and the rest under whatever directory it is started in, always under
        # those same names, so nine decks in one folder would overwrite each
        # other's results. The folder is also what the extractor looks for.
        write_deck(this_lance, terminal, tpa, name, mu, direction,
                   out / f"stage3_{name}" / f"stage3_{name}.k", opts)
        stroke_text = "-" if terminal is None else f"{terminal['stroke']:+.3f}"
        n_el = (len(this_lance["elements"])
                + (len(terminal["elements"]) if terminal else 0)
                + (len(tpa["elements"]) if tpa else 0))
        ramp = opts.get("t_ramp", T_RAMP_INSERT if direction < 0 else T_RAMP)
        print(f"{name:>20} {mu:>5.2f} "
              f"{'-' if gap is None else f'{gap:8.3f}':>8} "
              f"{stroke_text:>8} {opts.get('tooth_h', TOOTH_H):>8.3f} "
              f"{ramp * 1e3:>8.1f} {n_el:>9}")
        index["cases"][name] = {"mu": mu, "tpa_clearance": gap,
                                "direction": direction,
                                "stroke": terminal["stroke"] if terminal else None,
                                "t_ramp": ramp,
                                "tooth_h": opts.get("tooth_h", TOOTH_H),
                                "baseline": ("extract_mu020"
                                             if name.startswith("extract_mu020_")
                                             else None)}
        # A run that builds its own lance has its own node numbering, so it
        # needs its own map. One shared map was written for all nine the first
        # time round, and the fine-mesh run was then read at the wrong places
        # entirely - node 1766 is the lance tip in the default mesh and sits at
        # x = 6.96 in the refined one. Its forces were unaffected, because the
        # root nodes come first and keep their IDs, which is exactly what made
        # the error hard to see.
        if this_lance is not lance:
            index["cases"][name]["nodes"] = node_map(this_lance)

        if direction == 0:
            # The regression run. Prescribed tip lift, nothing else touching
            # the lance, so the root reaction is the Stage 1 tip force with the
            # tooth's stiffening in it. Carried here so the extractor compares
            # against a number written before the run, not after it.
            #
            # The baseline is Stage 1's own FE result at the same mesh, not the
            # beam-theory value. Stage 1 measured a 4.2 % model-form gap
            # between the 3D FE and Timoshenko and put it on record; that gap
            # belongs to the beam idealisation and is still there in Stage 3.
            # Checking the tooth against beam theory would charge this run for
            # a difference Stage 1 already accounted for.
            f_tooth = beam_flexibility(G, M, L, L, with_tooth=True)
            f_plain = beam_flexibility(G, M, L, L, with_tooth=False)
            stiffening = f_plain / f_tooth
            s1 = stage1_reference()
            index["cases"][name].update({
                "tip_lift": Y_TIP,
                "P_beam_theory": Y_TIP / f_tooth,
                "P_stage1_fe": s1,
                "P_predicted": (s1 * stiffening) if s1 else Y_TIP / f_tooth,
                "baseline_is": ("Stage 1 FE at the same mesh, times the tooth "
                                "stiffening" if s1 else "beam theory only - "
                                "Stage 1 results not found"),
                "P_stage1_beam": tip_force_timoshenko(G, M),
                "tooth_stiffening_pct": 100.0 * (stiffening - 1.0)})

    (REPO / "results").mkdir(exist_ok=True)
    (REPO / "results" / "stage3_model.json").write_text(json.dumps(index, indent=2))
    write_run_scripts(out, [c[0] for c in CASES])
    print(f"\nDecks in {out.relative_to(REPO)}, one folder each")
    print("Model index in results/stage3_model.json")
    print("run_all.bat / run_all.sh written next to the folders")


if __name__ == "__main__":
    main()
