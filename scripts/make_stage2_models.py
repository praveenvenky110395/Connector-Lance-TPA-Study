#!/usr/bin/env python3
"""
Writes the Stage 2 LS-DYNA decks.

Stage 2 adds exactly one thing to Stage 1: contact. The lance is the Stage 1
lance, imported from the Stage 1 generator so it cannot drift - same length,
width, thickness, material, mesh and root clamp. What is new is a rigid plate
above the tip with an inclined lower face, driven along -x. As it advances by s
the face descends at the tip by s tan(angle), so the plate pushes the lance down
without anything being prescribed on the lance itself.

Why the inclination is on the rigid plate and not on the lance: the contact sits
on the lance's tip edge, so the contact normal is the plate's normal. The plate
does not rotate, so the normal does not either. Put the wedge face on the lance
instead and the 6.4 deg tip rotation from Stage 1 would tilt the contact normal
with it - and at 30 deg that is not a small correction.

What gets compared:

    W/P = (tan(angle) + mu) / (1 - mu tan(angle))

W is the summed axial reaction at the root, P the summed transverse one, both
from spcforc at the same instant. Global equilibrium of the lance makes those
two the axial and transverse components of the contact force. The ratio holds no
stiffness at all, so the Stage 1 difference between the FE model and the beam
solution cancels out of it and Stage 2 tests only the contact and the friction.

Units: mm, tonne, s, N, MPa. Geometry is frozen in docs/numerical_spec.md.
"""

import math
from pathlib import Path

from make_stage1_models import (
    build_mesh, check_mesh, load_curve,
    L, B, T, Y_TIP, E, NU, RHO, NY,
    T_RAMP, T_END, T_CURVE,
)

REPO = Path(__file__).resolve().parents[1]

# Mesh level. n4 is the Stage 1 comparison mesh and Stage 2 stays on it - the
# question here is the contact, not the discretisation of the beam.
N_THRU = 4

# Rigid plate. The modulus here deforms nothing - the part is rigid - it only
# sizes the contact penalty stiffness, and that is the reason it is not steel.
#
# With a steel modulus against PBT the penalty comes out about 21 times stiffer
# than the parts in contact, and LS-DYNA then asks for a time step around 1.9e-8
# against the 5.3e-8 the lance actually sets. The penalty spring, not the
# structure, would be dictating the step. Setting the plate's modulus to the
# lance's puts both sides of the contact on the same stiffness and lifts the
# contact limit by sqrt(210000/9800), a factor of 4.6, which clears it.
PLATE_RHO = 7.85e-9
PLATE_E = E
PLATE_NU = 0.30

# Initial clearance between the plate face and the lance tip, mm. Small but not
# zero: starting in exact contact makes the first cycle report an initial
# penetration, and the run then opens with a force spike that has nothing to do
# with the mechanics. The travel below carries the gap, so the deflection at the
# end of the ramp is still exactly Y_TIP.
GAP = 0.010

# Plate mesh. It is rigid, so this only has to be fine enough for the contact
# search to find a segment under every slave node.
PLATE_DX = 0.10
PLATE_NZ = 2
PLATE_OVERHANG = 0.25   # each side in y, so the plate is wider than the lance
PLATE_BACK = 0.30       # how far the plate reaches behind the first contact
PLATE_FRONT = 0.30      # and ahead of the last

NODE_OFFSET = 2000      # plate node and element ids start above the lance's
ELEM_OFFSET = 2000

# The runs. 30 deg is the lead-in angle, 45 deg the return face; both are
# measured from the insertion axis. The friction sweep sits on the lead-in
# angle, and 45 deg gets the frictionless case because W/P = 1.000 exactly
# there - that is the regression test against the Stage 1 spring force.
CASES = [
    (30.0, 0.00), (30.0, 0.10), (30.0, 0.20), (30.0, 0.30),
    (45.0, 0.00), (45.0, 0.20),
]


def plate_travel(angle_deg):
    """Travel needed to close the gap and then deflect the lance by Y_TIP."""
    return (Y_TIP + GAP) / math.tan(math.radians(angle_deg))


def build_plate(angle_deg, first_node, first_elem):
    """
    Structured hex block whose lower face is the inclined plane

        z = A - x tan(angle),    A = T + GAP + L tan(angle)

    so that at x = L the face starts GAP above the lance tip and, after the
    plate has moved -x by s, sits at T + GAP - s tan(angle). The face therefore
    descends with x, which is what makes a plate moving towards the root push
    the tip down rather than lift it.
    """
    ta = math.tan(math.radians(angle_deg))
    travel = plate_travel(angle_deg)

    x0 = L - PLATE_BACK
    x1 = L + travel + PLATE_FRONT
    nx = max(2, round((x1 - x0) / PLATE_DX))
    ny = NY                                  # same divisions as the lance in y
    nz = PLATE_NZ

    a_const = T + GAP + L * ta
    def z_low(x):
        return a_const - x * ta

    # Flat top, set above the highest point of the inclined face.
    z_top = z_low(x0) + 0.40

    def nid(i, j, k):
        return first_node + i * (ny + 1) * (nz + 1) + j * (nz + 1) + k

    nodes, elements = [], []
    for i in range(nx + 1):
        x = x0 + (x1 - x0) * i / nx
        lo = z_low(x)
        for j in range(ny + 1):
            y = -PLATE_OVERHANG + (B + 2 * PLATE_OVERHANG) * j / ny
            for k in range(nz + 1):
                nodes.append((nid(i, j, k), x, y, lo + (z_top - lo) * k / nz))

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                elements.append((first_elem + len(elements), [
                    nid(i, j, k), nid(i + 1, j, k),
                    nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                ]))

    return {"nodes": nodes, "elements": elements, "travel": travel,
            "x0": x0, "x1": x1, "z_top": z_top, "face_at_tip": z_low(L)}


def write_deck(lance, plate, angle_deg, mu, path):
    """One LS-DYNA keyword deck. Fields are 10 characters wide."""
    ratio = (math.tan(math.radians(angle_deg)) + mu) / \
            (1.0 - mu * math.tan(math.radians(angle_deg)))

    lines = [
        f"$ Stage 2. Rigid plate inclined {angle_deg:.0f} deg from the insertion axis,",
        f"$ Coulomb friction mu = {mu:.2f}, plate travel {plate['travel']:.4f} mm.",
        f"$",
        f"$ Expected W/P = (tan a + mu)/(1 - mu tan a) = {ratio:.4f}",
        f"$ with W = summed x reaction and P = summed z reaction at the root,",
        f"$ both from spcforc at the same instant, taken during the hold.",
        f"$",
        f"$ Lance geometry and material are identical to Stage 1. Units mm-tonne-s-N-MPa.",
        "*KEYWORD",
        "*TITLE",
        f"Stage2 a{angle_deg:.0f} mu{mu:.2f} n{N_THRU}",

        "*CONTROL_TERMINATION",
        f"{T_END:10.4E}",

        "*CONTROL_TIMESTEP",
        "$#  dtinit    tssfac      isdo    tslimt     dt2ms",
        "       0.0       0.9         0       0.0       0.0",
        "$ dt2ms = 0, so no mass scaling. The rigid plate does not enter the",
        "$ time step calculation - the lance sets it, as in Stage 1.",

        "*CONTROL_ENERGY",
        "$#    hgen      rwen    slnten     rylen",
        "         2         2         2         2",
        "$ slnten must be 2 or glstat reports no sliding interface energy, which",
        "$ is the one number that tells whether the contact is injecting energy.",
    ]

    for name in ("GLSTAT", "MATSUM", "SPCFORC", "RCFORC", "NODOUT"):
        lines += [f"*DATABASE_{name}", f"{1.0e-5:10.4E}"]
    lines += ["*DATABASE_BINARY_D3PLOT", f"{2.5e-4:10.4E}"]

    lines += [
        "*PART",
        "Lance",
        "$#     pid     secid       mid",
        "         1         1         1",
        "*SECTION_SOLID",
        "$#   secid    elform",
        "         1        -1",
        "*MAT_ELASTIC",
        "$#     mid        ro         e        pr",
        f"         1{RHO:10.3E}{E:10.1f}{NU:10.2f}",

        "*PART",
        "Plate",
        "$#     pid     secid       mid",
        "         2         2         2",
        "*SECTION_SOLID",
        "$#   secid    elform",
        "         2         1",
        "*MAT_RIGID",
        "$#     mid        ro         e        pr",
        f"         2{PLATE_RHO:10.3E}{PLATE_E:10.1f}{PLATE_NU:10.2f}",
        "$ cmo = 1 with con1 = 5 locks y and z, con2 = 7 locks all rotations.",
        "$ x is left free because x is what the prescribed motion drives.",
        "$#     cmo      con1      con2",
        "       1.0         5         7",
        "$ Third card, and it is not optional. *MAT_RIGID expects it even when",
        "$ there is no local system to define; leave it out and the reader walks",
        "$ into the next keyword and stops with a format error on this card.",
        "$#     lco        a1        a2        a3        v1        v2        v3",
        "         0       0.0       0.0       0.0       0.0       0.0       0.0",

        "*NODE",
    ]
    for n, x, y, z in lance["nodes"]:
        lines.append(f"{n:8d}{x:16.8f}{y:16.8f}{z:16.8f}")
    for n, x, y, z in plate["nodes"]:
        lines.append(f"{n:8d}{x:16.8f}{y:16.8f}{z:16.8f}")

    lines.append("*ELEMENT_SOLID")
    for eid, conn in lance["elements"]:
        lines.append(f"{eid:8d}       1")
        lines.append("".join(f"{n:8d}" for n in conn))
    for eid, conn in plate["elements"]:
        lines.append(f"{eid:8d}       2")
        lines.append("".join(f"{n:8d}" for n in conn))

    for set_id, node_ids in ((1, lance["root"]), (3, lance["spine"])):
        lines += ["*SET_NODE_LIST", f"{set_id:10d}"]
        for i in range(0, len(node_ids), 8):
            lines.append("".join(f"{n:10d}" for n in node_ids[i:i + 8]))

    lines += [
        "*DATABASE_HISTORY_NODE_SET",
        "$ The node line from root to tip. Its tip node is the contact point, and",
        "$ its z displacement must reach -0.600 at the end of the ramp.",
        "         3",

        "*BOUNDARY_SPC_SET",
        "$ Root face, all degrees of freedom - unchanged from Stage 1.",
        "$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz",
        "         1         0         1         1         1         1         1         1",

        "*CONTACT_AUTOMATIC_ONE_WAY_SURFACE_TO_SURFACE",
        "$ Deformable lance as slave, rigid plate as master. One-way, so slave",
        "$ nodes are checked against master segments and not the reverse - which is",
        "$ what you want when the master is rigid.",
        "$",
        "$ Not a FORMING contact. Those are written for shell blanks and reject",
        "$ solid elements outright; this lance is meshed with hexes.",
        "$#    ssid      msid     sstyp     mstyp",
        "         1         2         3         3",
        "$ fs = fd and dc = 0, so friction is plain Coulomb with no velocity",
        "$ dependence - which is what the closed-form relation assumes. vdc = 0:",
        "$ no contact damping, because damping would add force to the quantity",
        "$ being measured.",
        "$#      fs        fd        dc        vc       vdc    penchk",
        f"{mu:10.2f}{mu:10.2f}       0.0       0.0       0.0         0",
        "$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf",
        "       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0",

        "*BOUNDARY_PRESCRIBED_MOTION_RIGID",
        "$ Plate driven along -x. The inclined face converts that into the",
        "$ transverse deflection of the lance - nothing is prescribed on the lance.",
        "$#     pid       dof       vad      lcid        sf",
        "         2         1         2         1       1.0",

        "*DEFINE_CURVE",
        "         1",
    ]
    for t, d in load_curve():
        # load_curve() ramps to -Y_TIP; rescale it to the plate travel in -x.
        lines.append(f"{t:20.10E}{d / Y_TIP * plate['travel']:20.10E}")

    lines.append("*END")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def check_plate(lance, plate, angle_deg):
    """
    Geometry checks that have to hold before the deck is worth running.

    The one that matters is clearance: the plate must touch the lance only at
    the tip edge and must not run into the top face further back. With the tip
    deflected by w(L) and the plate advanced to match, the gap at station x is

        (L - x) tan(angle) - w(L) + w(x)

    which stays positive as long as tan(angle) exceeds the tip slope, 0.1125 rad
    here. Worth checking rather than asserting, since it is the difference
    between a wedge and a plate ploughing through the beam.
    """
    ta = math.tan(math.radians(angle_deg))
    r = 1.5 * Y_TIP / L                      # tip slope for a tip-loaded beam

    if ta <= r:
        raise ValueError(f"angle {angle_deg} deg is shallower than the tip slope")

    # Face just above the tip at t = 0, by exactly the gap.
    if abs(plate["face_at_tip"] - (T + GAP)) > 1e-12:
        raise ValueError("plate face does not start one gap above the tip")

    # Clearance along the lance at full travel. The tip itself is excluded: the
    # clearance there is zero by construction, because that is the contact. What
    # is being tested is the top face behind it, so the scan stops one element
    # short of the tip and that element's clearance is what gets reported.
    def shape(xi):
        return (3.0 * xi**2 - xi**3) / 2.0

    def clearance(x):
        return (L - x) * ta - Y_TIP + Y_TIP * shape(x / L)

    back = L / max(1, round(L / (T / N_THRU)))       # one element length in x
    worst, worst_x = None, None
    for i in range(1, 201):
        x = (L - back) * i / 200.0
        g = clearance(x)
        if worst is None or g < worst:
            worst, worst_x = g, x
    if worst < 0.0:
        raise ValueError(f"plate would cut into the lance at x = {worst_x:.3f} mm")

    # Plate must still be clear of the root end after the full stroke.
    if plate["x0"] - plate["travel"] < 0.5 * L:
        raise ValueError("plate reaches too far back over the lance")

    return worst, worst_x


def main():
    lance = build_mesh(N_THRU)
    check_mesh(lance)
    # Stage 1 appends a control node for its rigid tip section. Stage 2 loads the
    # tip through contact instead, so that node is dropped and the lance mesh is
    # exactly the Stage 1 grid, node numbering included.
    lance["nodes"] = lance["nodes"][:-1]

    print(f"Lance: n{N_THRU}, {len(lance['nodes']):,} nodes, "
          f"{len(lance['elements']):,} elements  (Stage 1 mesh, unchanged)")

    written = []
    print(f"\n{'deck':>26} {'angle':>7} {'mu':>6} {'travel':>9} "
          f"{'W/P target':>11} {'clearance':>10}")
    for angle_deg, mu in CASES:
        plate = build_plate(angle_deg, len(lance["nodes"]) + NODE_OFFSET,
                            len(lance["elements"]) + ELEM_OFFSET)
        worst, worst_x = check_plate(lance, plate, angle_deg)

        ta = math.tan(math.radians(angle_deg))
        ratio = (ta + mu) / (1.0 - mu * ta)

        name = f"stage2_a{angle_deg:.0f}_mu{mu * 100:03.0f}.k"
        path = REPO / "ls-dyna/stage2" / name
        write_deck(lance, plate, angle_deg, mu, path)
        written.append(path)

        print(f"{name:>26} {angle_deg:>6.0f}° {mu:>6.2f} "
              f"{plate['travel']:>9.4f} {ratio:>11.4f} {worst:>10.4f}")

    print(f"\n{len(written)} decks written")
    for path in written:
        print(f"  {path.relative_to(REPO)}")

    print(f"""
Run order
  1  a45 mu 0.00     the regression test. W/P must come out 1.000.
  2  a30 mu 0.00     contact formulation at the second angle, still no friction
  3  a30 mu 0.10 / 0.20 / 0.30      the friction sweep
  4  a45 mu 0.20     friction at the second angle

If d3hsp still asks for a smaller time step than the lance sets, the levers in
order are: tssfac 0.9 -> 0.6 in *CONTROL_TIMESTEP, then sfm 1.0 -> 0.5 on the
third contact card. Both are reportable settings, not fixes to hide - whichever
one is used goes in the report with the reason.

Check before recording any force:
  - tip node z displacement at the end of the ramp. It should read -0.600, and
    slightly less is expected - penalty contact lets the tip penetrate a little.
    Record what it reads; the ratio is what the comparison rests on, and the
    ratio does not care.
  - sliding interface energy in glstat against internal energy. Contact that
    injects energy shows up here and nowhere else.
  - rcforc against spcforc. The contact resultant and the root reaction are the
    same force seen from two sides, so they must agree in x and in z.

Extract, from spcforc at the root, averaged over 5.25 to 5.50 ms:
  W = summed x reaction      P = summed z reaction      ratio = W/P

The ratio carries no stiffness, so the Stage 1 deviation from beam theory does
not propagate into it. If W/P is right and Stage 1 was 3.7 % off, both
statements are true at once and neither contradicts the other.
""")


if __name__ == "__main__":
    main()
