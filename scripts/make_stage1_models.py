#!/usr/bin/env python3
"""
Writes the Stage 1 LS-DYNA decks.

Stage 1 is a cantilever with a prescribed tip deflection. No contact, linear
elastic. It's meshed here instead of by hand because the four mesh levels are
then one loop, and because anyone can regenerate the exact same model.

The tip deflection goes on a rigid tip section through a control node on the
neutral axis, not on the individual face nodes. The reason is in the comment on
the *CONSTRAINED_NODAL_RIGID_BODY card and in docs/numerical_spec.md.

Units: mm, tonne, s, N, MPa. Geometry is frozen in docs/numerical_spec.md.
"""

from pathlib import Path

# Lance, mm. x runs along the length, y across the width, z through the
# thickness. Bending happens in z.
L = 8.00
B = 2.50
T = 0.80
Y_TIP = 0.60          # prescribed tip deflection, downwards

# BASF Ultradur B 4300 G6, PBT-GF30. E from ISO 527, density from ISO 1183.
E = 9800.0
NU = 0.35
RHO = 1.53e-9

# Ramp, termination and curve end are three separate times.
#
#   0 -------- T_RAMP -------- T_END -------- T_CURVE
#   |  smoothstep  |   hold    |   hold       |
#                              ^ analysis stops here
#
# The 5 ms ramp is about 26 first-bending periods (T1 = 0.196 ms). Quasi-static
# behaviour is checked from the KE/IE response in glstat.
#
# Termination sits inside the hold, and the curve runs past termination, so the
# prescribed displacement is defined on every step the solver takes including
# the last. The 0.5 ms hold also serves as a settling check: the reaction force
# should be flat across it.
T_RAMP = 5.0e-3       # smoothstep ramp completes
T_END = 5.5e-3        # analysis terminates, inside the hold
T_CURVE = 6.0e-3      # load curve continues past termination

# Elements across the width. Held constant - the load case has no gradient along
# y, so this does not refine with the through-thickness count. Not 1, because the
# section carries some anticlastic curvature from the Poisson effect. Set to 6
# rather than 3 for element aspect ratio at the finest mesh: at n8 the element is
# 0.1 mm through thickness, giving 8:1 across the width at NY=3 and 4:1 at NY=6.
NY = 6

# Stations along the length where the deflected shape is read, mm. Quarter
# points, so the comparison covers the whole beam rather than just the tip.
STATIONS = (2.0, 4.0, 6.0, 8.0)

# Point mass and rotational inertia carried by the tip control node, in tonne and
# tonne mm^2. The rigid tip section is a flat plane of nodes and its rotation
# about the bending axis is free, so the inertia it gets from the mesh alone is
# tiny - small enough that the explicit solver error terminates on it. These
# values regularise that. They are four orders below the beam mass of 2.45e-8 t,
# and they are held fixed across all four meshes so the only thing changing
# through the convergence sequence is the mesh. Being inertial terms, they affect
# the transient and not the settled force that gets recorded.
CONTROL_MASS = 1.0e-10
CONTROL_INERTIA = 1.0e-9
MASS_EID = 90001
INERTIA_EID = 90002

REPO = Path(__file__).resolve().parents[1]

# Formulations compared on a common mesh (n4). n4 is the comparison mesh,
# not a declared converged mesh - convergence is assessed in post-processing.
FORMULATIONS = {
    -1: "fully integrated S/R, handles poor aspect ratio - baseline",
     1: "one-point integration, needs hourglass control",
     2: "fully integrated, expect some stiffening in bending",
}


def build_mesh(n_thru):
    """Structured hex grid with n_thru elements through the thickness."""
    h = T / n_thru
    nx = max(1, round(L / h))
    nz = n_thru

    def node_id(i, j, k):
        return i * (NY + 1) * (nz + 1) + j * (nz + 1) + k + 1

    nodes = []
    for i in range(nx + 1):
        for j in range(NY + 1):
            for k in range(nz + 1):
                nodes.append((node_id(i, j, k),
                              L * i / nx, B * j / NY, T * k / nz))

    # Hex node order: 1-4 on the lower z face, 5-8 directly above. Getting this
    # backwards flips the Jacobian, which is what check_mesh catches - it caught
    # exactly that while this was being written.
    elements = []
    for i in range(nx):
        for j in range(NY):
            for k in range(nz):
                elements.append((len(elements) + 1, [
                    node_id(i, j, k),
                    node_id(i + 1, j, k),
                    node_id(i + 1, j + 1, k),
                    node_id(i, j + 1, k),
                    node_id(i, j, k + 1),
                    node_id(i + 1, j, k + 1),
                    node_id(i + 1, j + 1, k + 1),
                    node_id(i, j + 1, k + 1),
                ]))

    root = sorted(node_id(0, j, k) for j in range(NY + 1) for k in range(nz + 1))
    tip = sorted(node_id(nx, j, k) for j in range(NY + 1) for k in range(nz + 1))

    # A line of nodes running root to tip along the top surface at mid width.
    # This is what nodout writes, and it gives the deflected shape in one set at
    # every mesh level. Mid width because NY is even; top surface because the
    # neutral plane is only a node row when nz is even, and n1 has nz = 1.
    spine = [node_id(i, NY // 2, nz) for i in range(nx + 1)]

    # Control node for the tip, sitting on the neutral axis at the centre of the
    # tip face. No solid element touches it, so it adds no stiffness - it exists
    # to be the point the deflection is prescribed at and read from, and to carry
    # the point mass and inertia the rigid body needs. On the even meshes it lands
    # on top of an existing node; coincident nodes in the same rigid body are
    # harmless.
    control = len(nodes) + 1
    nodes.append((control, L, B / 2.0, T / 2.0))

    return {"h": h, "nz": nz, "nx": nx, "nodes": nodes, "elements": elements,
            "root": root, "tip": tip, "spine": spine, "control": control}


def check_mesh(mesh):
    """Every element must have positive volume, and they must sum to L*B*T."""
    coords = {n[0]: n[1:] for n in mesh["nodes"]}

    def tet_volume(a, b, c, d):
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        w = [d[i] - a[i] for i in range(3)]
        return (u[0] * (v[1] * w[2] - v[2] * w[1])
                - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0])) / 6.0

    # A hex splits into these five tetrahedra.
    tets = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6), (3, 4, 6, 7), (1, 4, 5, 6)]

    total = 0.0
    for _, connectivity in mesh["elements"]:
        p = [coords[n] for n in connectivity]
        volume = sum(tet_volume(p[a], p[b], p[c], p[d]) for a, b, c, d in tets)
        if volume <= 0:
            raise ValueError("negative Jacobian - check the node ordering")
        total += volume

    expected = L * B * T
    if abs(total - expected) > 1e-9 * expected:
        raise ValueError(f"mesh volume {total:.6f}, expected {expected:.6f}")


def load_curve(n_points=501):
    """
    Displacement against time: quintic smoothstep to T_RAMP, then held to T_CURVE.

    Both the first and second derivatives of the smoothstep are zero at each end.
    A straight ramp would start with a velocity step, which rings the cantilever -
    exactly the thing the quasi-static energy check is there to rule out.

    501 points rather than 101. The curve is interpolated linearly between points,
    so the prescribed velocity steps at every one of them and each step is a small
    kick. Five times as many points makes each kick five times smaller and puts
    the excitation well above anything the structure responds to, which shows up
    as a flatter reaction force during the hold.

    The hold runs past termination so the prescribed displacement is defined on
    every step the solver takes. An earlier version ended the curve exactly at
    termination, which left the last step undefined and distorted the tip elements
    in the final d3plot state while the rest of the model was unaffected.
    """
    points = []
    for i in range(n_points):
        s = i / (n_points - 1)
        shape = 10 * s**3 - 15 * s**4 + 6 * s**5
        points.append((s * T_RAMP, -Y_TIP * shape))
    points.append((T_CURVE, -Y_TIP))
    return points


def write_deck(mesh, elform, path):
    """Write one LS-DYNA keyword deck. Fields are 10 characters wide."""
    lines = [
        f"$ Stage 1 cantilever. {mesh['nz']} elements through thickness, "
        f"h = {mesh['h']:.3f} mm.",
        f"$ ELFORM {elform}: {FORMULATIONS[elform]}",
        "$ Linear elastic, no contact. Units mm-tonne-s-N-MPa.",
        "*KEYWORD",
        "*TITLE",
        f"Stage1 n{mesh['nz']} elform{elform}",

        "*CONTROL_TERMINATION",
        f"{T_END:10.4E}",

        "*CONTROL_TIMESTEP",
        "$#  dtinit    tssfac      isdo    tslimt     dt2ms",
        "       0.0       0.9         0       0.0       0.0",
        "$ dt2ms = 0, so no mass scaling. Check the cycle count in d3hsp before",
        "$ putting that claim in the report.",

        "*CONTROL_ENERGY",
        "$#    hgen      rwen    slnten     rylen",
        "         2         2         2         2",
        "$ hgen must be 2. On the default, glstat reports zero hourglass energy,",
        "$ which is not the same as there being none.",
    ]

    for name in ("GLSTAT", "MATSUM", "SPCFORC", "BNDOUT", "NODOUT"):
        lines += [f"*DATABASE_{name}", f"{1.0e-5:10.4E}"]
    lines += ["*DATABASE_BINARY_D3PLOT", f"{2.5e-4:10.4E}"]

    # The hourglass card only takes effect if the part points at it, so hgid on
    # the part card has to match. Leaving it blank silently falls back to the
    # default control and the *HOURGLASS card does nothing.
    hgid = 1 if elform == 1 else 0

    lines += [
        "*PART",
        "Lance",
        "$#     pid     secid       mid     eosid      hgid",
        f"         1         1         1         0{hgid:10d}",
        "*SECTION_SOLID",
        f"         1{elform:10d}",
    ]

    if elform == 1:
        # One-point integration has zero-energy modes, so it needs hourglass
        # control. Type 6 is the assumed-strain stiffness form, which is the one
        # meant for bending. qm = 1.0 rather than the 0.1 used with the viscous
        # forms: type 6 recovers the element bending stiffness, and scaling that
        # down by ten leaves the element too soft for the comparison this run is
        # supposed to make.
        lines += [
            "*HOURGLASS",
            "$#    hgid       ihq        qm",
            "         1         6       1.0",
        ]

    lines += [
        "*MAT_ELASTIC",
        "$#     mid        ro         e        pr",
        f"         1{RHO:10.3E}{E:10.1f}{NU:10.2f}",
        "*NODE",
    ]
    for n, x, y, z in mesh["nodes"]:
        lines.append(f"{n:8d}{x:16.8f}{y:16.8f}{z:16.8f}")

    lines.append("*ELEMENT_SOLID")
    for eid, connectivity in mesh["elements"]:
        lines.append(f"{eid:8d}       1")
        lines.append("".join(f"{n:8d}" for n in connectivity))

    tip_body = [mesh["control"]] + mesh["tip"]
    for set_id, node_ids in ((1, mesh["root"]), (2, mesh["tip"]),
                             (3, mesh["spine"]), (4, tip_body)):
        lines += ["*SET_NODE_LIST", f"{set_id:10d}"]
        for i in range(0, len(node_ids), 8):
            lines.append("".join(f"{n:10d}" for n in node_ids[i:i + 8]))

    lines += [
        "*DATABASE_HISTORY_NODE_SET",
        "$ *DATABASE_NODOUT on its own writes a header and no data - it needs to be",
        "$ told which nodes to trace. Set 3 is the line of nodes from root to tip",
        "$ used for the deflected-shape check.",
        "         3",

        "*DATABASE_HISTORY_NODE",
        "$ Control node as well. Its z displacement is the tip deflection the",
        "$ analytical solution is written for, and it must read -0.600.",
        f"{mesh['control']:10d}",
    ]

    lines += [
        "*BOUNDARY_SPC_SET",
        "$ Root face, all degrees of freedom. This is the perfect clamp.",
        "$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz",
        "         1         0         1         1         1         1         1         1",

        "*CONSTRAINED_NODAL_RIGID_BODY",
        "$ Tip face plus the control node, tied into one rigid body.",
        "$",
        "$ Beam theory says plane sections remain plane, and this is that sentence",
        "$ written as a constraint. It matters because the tip rotates 6.4 deg here:",
        "$ a plane section rotated by theta projects shorter onto z by t(1-cos theta),",
        "$ so prescribing the same z displacement at every node of the face - the",
        "$ obvious way to do it - forbids that projection change and forces a",
        "$ through-thickness strain of theta^2/2 into the tip elements. At this",
        "$ deflection that is 0.63 %, against a root bending strain of 1.125 %.",
        "$",
        "$ The rigid body has no rotational constraint, so the section is free to",
        "$ rotate. Only the z translation is prescribed.",
        "$#     pid       cid      nsid     pnode      iprt    drflag    rrflag",
        "         2         0         4         0         0         0         0",

        "*ELEMENT_MASS",
        "$ Added to the rigid body, not replacing what it gets from the mesh.",
        "$#     eid        id            mass",
        f"{MASS_EID:8d}{mesh['control']:8d}{CONTROL_MASS:16.4E}",

        "*ELEMENT_INERTIA",
        "$ Without this the run error terminates: a flat plane of nodes free to",
        "$ rotate has almost no inertia about the bending axis, and the angular",
        "$ accelerations that follow are not integrable at this time step.",
        "$#     eid        id      csid",
        f"{INERTIA_EID:8d}{mesh['control']:8d}",
        "$#     ixx       ixy       ixz       iyy       iyz       izz",
        f"{CONTROL_INERTIA:10.3E}{0.0:10.1f}{0.0:10.1f}"
        f"{CONTROL_INERTIA:10.3E}{0.0:10.1f}{CONTROL_INERTIA:10.3E}",

        "*BOUNDARY_PRESCRIBED_MOTION_RIGID",
        "$ Applied to the rigid body, not to a node of it. The body's z translation",
        "$ follows the curve; x, y and all three rotations stay free.",
        "$#     pid       dof       vad      lcid        sf",
        "         2         3         2         1       1.0",

        "*DEFINE_CURVE",
        "         1",
    ]
    for t, d in load_curve():
        lines.append(f"{t:20.10E}{d:20.10E}")

    lines.append("*END")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main():
    meshes = {}
    print(f"{'mesh':>5} {'h [mm]':>8} {'nodes':>8} {'elements':>9} {'control node':>13}")
    for n in (1, 2, 4, 8):
        mesh = build_mesh(n)
        check_mesh(mesh)
        meshes[n] = mesh
        print(f"{'n' + str(n):>5} {mesh['h']:>8.3f} "
              f"{len(mesh['nodes']):>8,} {len(mesh['elements']):>9,} "
              f"{mesh['control']:>13,}")

    # Node IDs to pull out of nodout for the deflected-shape check. A station is
    # only listed where it lands exactly on a node: n1 has 0.8 mm spacing, so
    # 2 mm and 6 mm fall between nodes and there is nothing to read there.
    print("\nDeflected-shape stations, node IDs in set 3 (written to nodout)")
    print(f"{'mesh':>5} " + "".join(f"{'x=' + f'{x:.0f}' + ' mm':>10}" for x in STATIONS))
    for n, mesh in meshes.items():
        cells = []
        for x in STATIONS:
            index = x / L * mesh["nx"]
            on_node = abs(index - round(index)) < 1e-9
            cells.append(f"{mesh['spine'][round(index)]:>10}" if on_node else f"{'-':>10}")
        print(f"{'n' + str(n):>5} " + "".join(cells))

    written = []

    # Mesh convergence: baseline formulation, all four levels.
    for n, mesh in meshes.items():
        path = REPO / "ls-dyna/stage1/mesh" / f"stage1_n{n}_elform-1.k"
        write_deck(mesh, -1, path)
        written.append(path)

    # Formulation comparison on the common mesh only. Running every formulation
    # at every mesh level would be 12 runs telling us what 6 already do.
    for elform in (1, 2):
        path = REPO / "ls-dyna/stage1/formulation" / f"stage1_n4_elform{elform}.k"
        write_deck(meshes[4], elform, path)
        written.append(path)

    print(f"\n{len(written)} decks written")
    for path in written:
        print(f"  {path.relative_to(REPO)}")

    print("""
Run order
  1  n1, n2, n4, n8 with ELFORM -1     mesh convergence
  2  n4 with ELFORM 1 and ELFORM 2     formulation comparison, common mesh

Ramp completes at 5.0 ms, analysis terminates at 5.5 ms. Read the force during
the hold, and check that it is flat across it - drift there means the response
has not settled.

Target is 3.6455 N (Timoshenko). Euler-Bernoulli is 3.6750 N; the 0.81 %
difference is transverse shear and is expected. Accept 3.573 to 3.718 N.

The tip force is the summed z reaction, read from SPCFORC at the root. Average
it over the hold rather than reading a single sample - the model is undamped and
the trace carries some ringing.

Check the control node in NODOUT first, before any force is recorded. Its z
displacement must read -0.600, which is what the analytical solution is written
for. Nodes off the neutral axis will read slightly more: the section has rotated
6.4 deg and a rotated plane projects shorter onto z. That is geometry, not error.

NODOUT also carries the deflected shape along the node line above, which is a
check on the response that does not involve the reaction force at all.

Expected behaviour:
  n1               outside the asymptotic range. A single linear element through
                   the thickness carries constant transverse shear and cannot
                   represent the parabolic distribution, and its aspect ratio is
                   the worst of the four. Excluded from the convergence fit.
  ELFORM -1        reference formulation for the mesh sequence
  ELFORM 1         check hourglass energy in glstat against internal energy
  ELFORM 2         stiffer response expected from full integration in bending

Record the measured values, including deviations from the expected behaviour.
""")


if __name__ == "__main__":
    main()
