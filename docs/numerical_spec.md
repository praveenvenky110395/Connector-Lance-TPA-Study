# Numerical Specification

**Connector Lance & TPA Study** — frozen inputs and analytical expectations.

| | |
|---|---|
| Version | 1.0 |
| Status | Inputs frozen. All three stages evaluated. Stage 3 force results quoted for five of eight runs; see Section 6 |
| Unit system | mm – tonne – s – N – MPa |
| Computed by | `scripts/analytical.py` (nothing in this document is hand-entered) |

> Every number in this document is produced by `scripts/analytical.py` from the inputs in
> its `Geometry` and `Material` blocks. The script re-derives each quantity by a second
> independent route and raises if the two disagree. Re-run it to regenerate these values.

---

## 1. Unit system

| Quantity | Unit | Note |
|---|---|---|
| Length | mm | |
| Mass | tonne (Mg) | |
| Time | s | |
| Force | N | |
| Stress | MPa | = N/mm² |
| Density | tonne/mm³ | 1 g/cm³ = 1e-9 t/mm³ |

---

## 2. Geometry

| Symbol | Value | Description |
|---|---|---|
| L | **8.00 mm** | Lance free length, root clamp to tooth centreline |
| b | **2.50 mm** | Lance width |
| t | **0.80 mm** | Lance thickness |
| y | **0.60 mm** | Required deflection |
| L/t | **10.0** | Slenderness — screening criterion for beam-theory applicability |
| I | **0.106667 mm⁴** | Second moment of area, b·t³/12 |
| α | **30°** | Lead-in ramp, insertion side |
| β | **45°** | Return face, retention side |

### 2.1 Angle convention — read this before opening CAD

**Both angles are measured from the insertion axis.** A shallow ramp is a *small* angle.

- α → 0° is a face parallel to the insertion direction: pure sliding, W → µP
- α → 90° is a face normal to the insertion direction: self-locking, cannot be inserted

Drawing these from the perpendicular instead inverts every Stage 2 number and Stage 1
will not reconcile with Stage 2. This is the single easiest way to lose a day.

### 2.2 Geometric consistency constraint

**Tooth protrusion = notch depth = y = 0.60 mm.** All three must be equal.

If the tooth protrudes less than the notch is deep, the lance never reaches the Stage 1
deflection, the µ = 0 regression test fails for a purely geometric reason, and the two
stages are not comparable.

### 2.3 Scale

These are **representative dimensions in the range of automotive connector locking
features**. No claim is made that they correspond to any particular terminal size class,
and the geometry is synthetic — it is not a reproduction of any manufacturer's product.

---

## 3. Material

**BASF Ultradur® B 4300 G6 (PBT-GF30)** — published datasheet values.

| Property | Value | Standard |
|---|---|---|
| Tensile modulus E | **9,800 MPa** | ISO 527-1/-2 |
| Stress at break | **137 MPa** | ISO 527-1/-2 |
| Strain at break | **3.0 %** | ISO 527-1/-2 |
| Density | **1,530 kg/m³** → 1.53e-9 t/mm³ | ISO 1183 |
| Tensile creep modulus, 1000 h / 23 °C | **7,500 MPa** | ISO 899-1 |
| Poisson's ratio ν | **0.35 — assumed** | not published |

### 3.1 Why PBT and not PA66

PA66's modulus changes substantially between dry and conditioned states. A part whose
function is holding a terminal at a position tolerance cannot have its stiffness depend
on ambient humidity. PBT's moisture uptake is far lower, which is why it dominates
connector housings. This is the reason for the choice, not "it is typical".

### 3.2 Poisson's ratio

ν is not on the datasheet and is assumed. Sweeping 0.30 → 0.40 moves the Timoshenko
force from 3.6466 N to 3.6444 N — a spread below 0.1 %. **The assumption is immaterial
to the result**, and this sensitivity is reported rather than the assumption hidden.

### 3.3 Permissible strain

| | |
|---|---|
| Strain at break (ISO 527) | 3.00 % |
| Factor, one-time assembly, glass-reinforced | **0.50** |
| **Permissible strain** | **1.50 %** |
| Design strain | 1.125 % |
| Utilisation | **75 %** |

0.50 is the conservative end of the published design-guide range for glass-reinforced
grades in one-time assembly.

> **OPEN ITEM.** Cite the specific design guide consulted. A referenced factor beats a
> plausible one, and the interviewer writes specifications for a living. At 0.60 the
> design still passes, at 62 % utilisation.

### 3.4 Creep and relaxation — what may and may not be claimed

**May state:** the datasheet reports a reduction from 9,800 MPa tensile modulus to
7,500 MPa tensile creep modulus after 1000 h at 23 °C (ISO 899-1).

**Must state:** the lance is held at approximately constant *strain*, so **stress
relaxation** — not creep — is the strictly correct quantity. The creep modulus is an
indicator that long-term stiffness degradation can be significant; it is not a
relaxation prediction.

**Must not state:** any retention force after 1000 h. No such number is claimed.

---

## 4. Analytical targets

### 4.1 Stage 1 — perfect root clamp, linear elastic, no contact

| Quantity | Value |
|---|---|
| Euler–Bernoulli, P_EB | **3.6750 N** |
| Shear correction | **0.8100 %** |
| **Timoshenko, P_T** | **3.6455 N** ← *analytical reference for the FE comparison* |
| Acceptance band (±2 %) | 3.573 – 3.718 N |
| Root stress (linear, from P_EB) | 110.25 MPa |
| Root strain | 1.1250 % |
| First bending mode f₁ | 5,110 Hz |
| T₁ | 0.1957 ms |
| Wave speed | 2,530,855 mm/s |

**κ convention: classic, κ = 5/6** (rectangular section). Cowper's
κ = 10(1+ν)/(12+11ν) gives 0.7925 % instead of 0.8100 % — a 0.02 % difference in the
final force. The choice is stated because the question will be asked, not because it
matters.

**P_T, not P_EB, is the analytical reference for the FE comparison.** The FE result is
expected to approach P_T with mesh refinement under the stated idealized assumptions. The
direction and magnitude of the expected deviation are known in advance: approximately
0.8 % below Euler–Bernoulli, from transverse shear.

#### Deflected shape — a second check on the same model

The tip deflection is *prescribed* in the FE model, so the shape of the beam between root
and tip is not implied by it. Comparing that shape against the closed form tests the
response without the reaction force entering the comparison at all — it is a second,
independent check on the same run, not a restatement of the first.

w(x)/w(L) = [(3ξ² − ξ³) + 2rξ] / [2(1 + r)],  ξ = x/L, r = shear correction ratio

| x [mm] | ξ | w Euler–Bernoulli [mm] | w Timoshenko [mm] | difference |
|---|---|---|---|---|
| 2.0 | 0.25 | 0.0516 | 0.0524 | 1.53 % |
| 4.0 | 0.50 | 0.1875 | 0.1884 | 0.48 % |
| 6.0 | 0.75 | 0.3797 | 0.3803 | 0.15 % |
| 8.0 | 1.00 | 0.6000 | 0.6000 | prescribed |

The two theories differ by about 1.5 % at quarter span, which is close to the resolution
of the mesh sequence. **This check discriminates beam behaviour from a wrong model; it
does not discriminate Timoshenko from Euler–Bernoulli.** Stated so the claim is not
overreached.

Read from `nodout`, node set 3 — a line of nodes from root to tip along the top surface at
mid width, written by every Stage 1 deck. The tip node of that set must reach exactly
−0.600 mm, which confirms the prescribed motion was delivered before any result is read.

### 4.2 Stage 2 — wedge relation

W = P (tan θ + µ) / (1 − µ tan θ)

| µ | W insertion (α = 30°) | W retention (β = 45°) |
|---|---|---|
| **0.00** | **2.1218 N** | **3.6750 N** |
| 0.10 | 2.6418 N | 4.4917 N |
| 0.20 | 3.2297 N | 5.5125 N |
| 0.30 | 3.8997 N | 6.8250 N |

**Run µ = 0 first.** The friction term vanishes and W = P·tan α is pure geometry, so a
mismatch can only come from the contact definition — one unknown instead of two. And
with β = 45°, tan β = 1, so the retention force returns **exactly P_EB**: a regression
test on Stage 1, hidden inside Stage 2, at no cost.

Friction is a **parametric assumption with a documented sensitivity study**, never a
measured material property. Results are reported as a band.

---

## 5. Stage definitions

### Stage 1 — Structural verification

**Tools: LS-DYNA and Python only.** A second FE code was considered and dropped. The
closed-form solution is already the independent reference — it is not a solver — and the
days a second code would have cost are worth more spent on Stage 2 contact.

**Meshing.** Generated by `scripts/make_stage1_models.py`, so the four mesh levels are one
loop and anyone can regenerate the exact model. The generator checks every element
Jacobian and the total volume. Stage 2/3 geometry is meshed properly once it has ramps and
contact surfaces.

Six elements across the width at every level. Held constant — the load case has no
gradient along y, so it does not refine with the through-thickness count. Not 1, because
the section carries some anticlastic curvature from the Poisson effect. Six rather than
three is driven by element aspect ratio at the finest mesh: at n8 the element is 0.1 mm
through thickness, which gives 8:1 across the width at three elements and 4:1 at six.

At b/t = 3.12 the section is narrow enough that beam theory with E is the appropriate
reference; a much wider section would stiffen toward E/(1−ν²), about 14 % at ν = 0.35.

| Mesh | h [mm] | Nodes | Elements |
|---|---|---|---|
| n1 | 0.800 | 155 | 60 |
| n2 | 0.400 | 442 | 240 |
| n4 | 0.200 | 1,436 | 960 |
| n8 | 0.100 | 5,104 | 3,840 |

The node counts include one control node at the tip, described below. No solid element
touches it, so it adds no stiffness.

Perfect root clamp (SPC all DOF on the root face), linear elastic, **no contact**, no
housing. Prescribed tip displacement y = 0.60 mm, quintic smooth-step, `DT2MS = 0`.

#### How the tip deflection is applied — and why not the obvious way

The tip face is tied into a **nodal rigid body** together with a control node on the
neutral axis at the centre of the face. The deflection is prescribed on the rigid body in
z only; x, y and all three rotations stay free.

The obvious alternative — prescribing the same z displacement at every node of the tip
face — was used first and is wrong on this geometry. The reason is geometric.

At y = 0.60 mm on an 8 mm cantilever the tip section rotates θ = 1.5 y / L = 0.1125 rad,
or **6.4°**. A plane section rotated by θ projects shorter onto the z axis by
t(1 − cos θ). Prescribing one z value across the whole face forbids that projection
change, so the section is forced to stretch through its thickness instead:

**ε_zz = θ² / 2 = 0.633 %**

against a root bending strain of 1.125 %. The constraint injects an artificial
through-thickness strain worth more than half the peak design strain of the model, at the
one location where the result is measured. Measured cost on the n4 mesh: about 7 %.

The rigid tip section removes it. Euler–Bernoulli assumes plane sections remain plane, and
a rigid end section is that assumption written as a constraint — the FE model is made to
adopt the kinematics the closed form is built on, which is what a verification study is
supposed to do. The cost is that the section cannot warp under transverse shear, which
Timoshenko theory allows; this is a small stiffening error in the opposite direction and
is accepted rather than removed, because the alternatives trade it for a less defensible
assumption.

The effect scales as θ², so it vanishes for small deflections. At δ/L = 0.075 it does not.

**Inertia on the control node.** The rigid tip section is a flat plane of nodes whose
rotation about the bending axis is free, and the inertia it inherits from the mesh about
that axis is very small. The explicit solver error terminates on it. A point mass of
1×10⁻¹⁰ t and a diagonal inertia of 1×10⁻⁹ t·mm² are added at the control node to
regularise it — four orders below the beam mass of 2.45×10⁻⁸ t, added to what the rigid
body gets from the mesh rather than replacing it, and held identical across all four
meshes so that the mesh is the only thing changing through the convergence sequence.
These are inertial terms: they affect the transient, not the settled force that is
recorded during the hold.

**Verification of the boundary condition itself.** The control node z displacement in
`nodout` must read −0.600 at termination. Nodes off the neutral axis read slightly more:
a node at Δz from the axis carries an additional Δz(cos θ − 1), which for the top surface
at this deflection is −0.0025 mm. That is the rotation projecting, not an error, and it is
checked rather than assumed.

#### Hypotheses tested and rejected

| Hypothesis | Test | Result |
|---|---|---|
| Shear locking in the fully integrated element | ELFORM −1, 1, 2 on a common mesh | 0.55 % spread across all three — not the element formulation |
| Over-constrained root clamp | lateral DOF released on the root face, one node held to prevent drift | no measurable change — not the root |

**Ramp, termination and curve end are three separate times.** The smooth-step ramp
completes at **5.0 ms** (≈ 26 × T₁, v ≈ 120 mm/s); the analysis terminates at **5.5 ms**,
inside the hold; the load curve continues to **6.0 ms**, past termination. The prescribed
displacement is therefore defined on every step the solver takes, including the last.

The 0.5 ms hold is also a **settling check**: the reaction force should be flat across it.
Drift over the hold means the response has not settled, which is a stronger statement than
the KE/IE ratio alone. The force is averaged over each half of the hold — 5.00–5.25 ms and
5.25–5.50 ms — rather than sampled at a single instant, because the model is undamped and
the reaction trace carries ringing. The difference between the two halves is the drift.

- Analytical reference: `scripts/analytical.py`
- LS-DYNA explicit, three formulations on a common mesh (n4):
  `ELFORM -1` (fully integrated S/R, poor-aspect-ratio formulation) as baseline;
  `ELFORM 1` (one-point integration with hourglass control, `IHQ=6, QM=1.0` — the
  assumed-strain stiffness form, at the coefficient that recovers the element bending
  stiffness rather than the 0.1 used with the viscous forms) — hourglass energy must
  stay below 10 % of internal;
  `ELFORM 2` (fully integrated) — expected to stiffen, since full integration
  struggles in bending.
- Mesh convergence on the baseline formulation only; formulation comparison on the
  common mesh only. Crossing the two would be 12 runs saying what 6 already say.
  n4 is the comparison mesh, not a declared converged mesh — convergence is assessed
  in post-processing from the observed order and the Richardson extrapolation.

**What the convergence study is for.** Not to show that a coarse mesh is inaccurate —
that is assumed, not demonstrated. Its purpose is to *justify the mesh carried into
Stages 2 and 3*, and to put a number on the discretization error.

Linear hex elements are nominally second-order accurate in displacement, so the error
is expected to fall roughly as 1/n². That is a prediction, made before running, in the
same way the Timoshenko offset is. From the three finest meshes the script computes:

- the **observed order of convergence** p = ln[(f₃−f₂)/(f₂−f₁)] / ln r, compared against
  the expected value of 2;
- the **Richardson extrapolation** to zero mesh size, f ≈ f₁ + (f₁−f₂)/(rᵖ−1), compared
  against the analytical solution;
- the **Grid Convergence Index** (ASME V&V 20 style, safety factor 1.25), which is an
  estimate of the discretization error on the finest mesh.

The GCI estimates discretization error; the difference from the analytical solution is a
separate quantity. They answer different questions and are reported separately.

The coarsest mesh is **excluded from the extrapolation** — one element through the
thickness is nowhere near the asymptotic range and would corrupt the fit. It is kept as
the coarse anchor of the convergence curve and as a smoke test: at 30 elements it runs
in seconds and confirms the deck is valid before the larger runs are committed to.
- Equilibrium: `SPCFORC` at the root must equal `BNDOUT` at the tip.

**No housing block in Stage 1.** Its local compliance would soften the response and
produce a discrepancy that is real physics rather than error.

**Do not cross the matrix:** mesh convergence on one formulation, formulation comparison
at one common mesh. Eight runs, not twenty.

### Stage 2 — Contact mechanics

**The plan was cut down, deliberately.** An earlier version of this section called for a
deformable terminal and a housing, `*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE` with `SOFT=1`
and `SOFT=2` compared, and contact damping of `VDC = 10–20`. That model would have added
contact, a second material, a second geometry and a softening housing all at once — and
the stage exists to verify **one** new effect. It was reduced to the smallest model that
still tests contact and friction against a closed form. The terminal and housing belong to
Stage 3, where there is no closed form to lose.

**Geometry.** The Stage 1 lance, unchanged — same length, width, thickness, material, mesh
(n4) and root clamp, with identical node numbering. Generated by importing the Stage 1
generator, so it cannot drift. What is added is a rigid plate above the tip whose lower
face is the inclined plane

z = T + gap + (L − x)·tan α

driven along −x. Advancing by s lowers that face at the tip by s·tan α, so the plate
deflects the lance without anything being prescribed on the lance at all. Travel to reach
y = 0.60 mm is (y + gap)/tan α — 1.0566 mm at 30°, 0.6100 mm at 45°, with a 0.010 mm
initial clearance so the run does not open on an initial penetration.

**Why the inclination is on the rigid body.** Contact sits on the lance's tip edge, so the
contact normal is the plate's normal, and the plate does not rotate. Put the wedge face on
the lance instead and the 6.4° tip rotation from Stage 1 would tilt the contact normal with
it; at α = 30° that is a 28 % error in tan α, not a correction. The same rotation that broke
the Stage 1 boundary condition would have broken Stage 2's reference.

**Clearance.** The plate must touch only the tip edge. At full travel the gap at station x is
(L − x)·tan α − y + y·w(x)/w(L), which stays positive as long as tan α exceeds the tip slope
1.5 y/L = 0.1125. The generator computes it and refuses to write a deck that fails: 0.093 mm
at 30°, 0.178 mm at 45°, one element behind the tip.

**Solver settings.**

- `*CONTACT_AUTOMATIC_ONE_WAY_SURFACE_TO_SURFACE`, deformable lance as slave, rigid plate
  as master. One-way, so slave nodes are checked against master segments and not the
  reverse — which is what is wanted when the master is rigid.
  **Not a `FORMING` contact.** Those read as the obvious choice for rigid tooling against
  a deformable part, and they are — for *shells*. They reject solid elements outright
  (`Error 40558`), and this lance is meshed with hexes.
- **Rigid plate modulus set to 9,800 MPa, the lance's own, not a steel value.** The plate
  is rigid, so its modulus deforms nothing; it only sizes the contact penalty stiffness.
  At a steel modulus the penalty comes out roughly 21× stiffer than the parts in contact,
  and LS-DYNA then requires a time step of about 1.9e-8 s against the 5.3e-8 s the lance
  itself sets — the penalty spring, not the structure, would be setting the step. Matching
  the moduli raises the contact limit by √(210000/9800) ≈ 4.6, to about 8.7e-8 s, which
  clears the structural step with a 1.6× margin.
- `FS = FD = µ`, `DC = 0` — plain Coulomb, no velocity dependence, which is what the closed
  form assumes.
- **`VDC = 0`. No contact damping.** Damping would add force to the quantity being measured.
  If a run needs damping to stay stable, that is a result to report, not a setting to bury.
- Default penalty scaling (`SFS = SFM = 1.0`).
- Rigid plate constrained through `*MAT_RIGID` `CMO = 1`, `CON1 = 5`, `CON2 = 7`: y, z and
  all rotations locked, x left free for the prescribed motion.
- Timing as Stage 1 — ramp to 5.0 ms, terminate 5.5 ms, curve to 6.0 ms.

**Runs.** α = 30° at µ = 0 / 0.10 / 0.20 / 0.30, and α = 45° at µ = 0 / 0.20. Six runs.

#### What Stage 2 is checked on, and why it is a ratio

W/P = (tan α + µ) / (1 − µ tan α)

with W the summed axial reaction at the root and P the summed transverse one, both from
`spcforc` at the same instant. Global equilibrium of the lance makes those two the axial and
transverse components of the contact force, so the relation follows from the contact normal
direction and the friction law alone.

**The ratio contains no stiffness** — not E, not I, not the deflection. So the Stage 1
difference between the FE model and the beam solution, which is a stiffness difference,
cancels out of it exactly. Stage 2 therefore tests only what Stage 2 adds. Stage 1 can sit
3.7 % above beam theory and Stage 2 can still be exact; both statements hold at once and
neither contradicts the other. This is what the staged structure was for.

**Regression test: α = 45°, µ = 0 gives W/P = 1.000 exactly.** No friction, no stiffness, no
material enters it — only the contact normal. Run this one first. If it is wrong the contact
is wrong and no other number in the stage is safe. Measured: **1.0002, +0.02 %.**

**Measurement window: 2.50–3.50 ms, while the plate is moving — not the hold.** This was
changed after the first full set of runs, and the reasoning is recorded because the change
itself is a result.

W/P = (tan α + µ)/(1 − µ tan α) is a *sliding* relation: it follows from Coulomb friction
with relative motion at the interface. The hold is the obvious window — the force is flat
there — but during the hold the plate is stationary, the contact is stuck, and a stuck penalty
interface carries whatever tangential force its stick spring held when motion stopped, less
whatever it sheds as the lance settles. Read in the hold the six runs give:

| µ | 0 (30°) | 0.10 | 0.20 | 0.30 | 0 (45°) | 0.20 (45°) |
|---|---|---|---|---|---|---|
| deviation | +0.03 % | −0.89 % | −1.52 % | −1.90 % | −0.07 % | −0.57 % |
| in standard errors | 0.1 | −4.0 | −4.8 | −4.8 | −0.1 | −3.8 |

Every friction case significant at about four standard errors, the miss growing with µ, and
both frictionless cases untouched — the signature of a tangential force relaxing, since with
µ = 0 there is none to relax. The frictional-work check had already said friction was fully
mobilised during sliding, coming out *above* the closed form rather than below.

The falsification condition was set before the data was examined: read while the plate moves,
the friction cases should land on the closed form and the µ-trend should vanish; if they
stayed ≈2 % low in every window, the friction model was at fault and the hold number would
stand. Result: worst deviation across all six drops to +0.41 %, with no µ-trend. The window
is not delicate — over all 35 windows between 0.8 and 1.5 ms long inside 2.2–4.2 ms, every
case has a median deviation within 0.2 % and no window in any case exceeds 0.93 %. 2.50–3.50
ms is simply the fastest part of the quintic stroke with the force already past half its
final value. Both readings are reported.

**Also extracted:** sliding interface energy against internal energy (`SLNTEN = 2`);
`RCFORC` against `SPCFORC`, which is the same force seen from two sides and must agree in x
and z; and the tip deflection reached. Two of those needed correcting once the first run was
read, and both corrections are to the criterion, not to the model:

- The 5 % screen on sliding interface energy only holds for the frictionless runs. There the
  counter contains nothing but penalty energy — the α = 45°, µ = 0 run gives 1.90 % — so it
  does measure how hard the contact works against itself. With friction on, the same counter
  also collects the frictional dissipation, which is real work: at α = 30°, µ = 0.3 the
  closed form puts it near 80 % of the internal energy. Screening that at 5 % would fail a
  run for behaving correctly. The friction cases are compared against
  E = ½ µ N s / cos α instead, with N = P / (cos α − µ sin α) and s = y / tan α — a second,
  independent check on the friction law, in energy rather than in force. The check that does
  hold for every case is total energy / initial energy, screened at ±1 %.
- KE/IE is read over the second half of the ramp only. Before that the lance holds almost no
  internal energy — 2.7 × 10⁻³ N·mm at 1 ms against 1.15 at the end — and the ratio is set
  by its own denominator rather than by anything dynamic. The same run reads 7.3 % at 1.02 ms
  and 0.32 % at 2.51 ms while the kinetic energy itself is still rising. The absolute peak
  and its time are printed for every run so the early part stays visible.

The tip deflection is close to 0.600 mm but not equal to it, and it can land either side —
+0.97 % on the α = 45°, µ = 0 run. The plate face is inclined, so the height of the face
under the tip depends on where the tip corner has moved to, and that corner both draws in
(the Stage 1 moment arm, 8.000 → 7.973 mm) and rides forward with the section rotation. This
is geometry, not penetration, and it is recorded rather than corrected — the ratio
does not depend on it.

### Stage 3 — the locking cycle

**Rewritten twice, and both rewrites are recorded because the reasoning is the content.**

The original plan asked for a TPA-on and TPA-off pull-out at 300 mm/s with no reference
beyond a limiting-case regression onto Stage 2. That does not work as written: with the
lance blocked, a force-driven pull-out in a linear-elastic model has nothing to limit it,
so the answer would be set by the contact penalty stiffness rather than by mechanics.

The second attempt modelled the TPA as a rigid prop under a prismatic lance and verified it
against the propped-cantilever closed form. Rigorous, and it produced a real result — a
prop part-way along the lance only stiffens it, it does not block release, and at 3.0 mm
with 0.10 mm clearance the retention force rises by just 17 %. But it is a beam exercise,
not a connector, and it left the part out.

**What is modelled now.** The lance carries an integral locking tooth: 45° retention face
on the root side, flat crest, 30° lead-in towards the tip, protruding 0.60 mm. The
protrusion is the lift needed to release — that is the design relation. A rigid terminal
runs underneath with a matching 45° shoulder and a recess the tooth sits in; a rigid TPA
sits above with a clearance. Both rigid parts take the lance's modulus, for the reason
Stage 2 established.

**The reference is Stage 1 × Stage 2**, with three named corrections:

| correction | size | why |
|---|---|---|
| tooth stiffens the lance | +1.1 % | thicker section, but where the moment is small |
| load on the tooth, not the tip | lever arm 6.26 mm not 8.00 mm | release condition set at the crest |
| faces rotate with the lance | **7.3–8.0°** | the faces are on the lance now |

The third is the finding. Lead-in 30° → 37.8° effective (steeper); retention 45° → 37.7°
(shallower). They converge. The insertion-to-retention asymmetry falls from about 2.2 on
the drawn angles to **1.30**, and it is 1.30 at every friction coefficient, because the
correction is geometric rather than frictional. **Design recommendation: draw the retention
face at 52.3° to end up with an effective 45°.**

Each force is reported as a bracket, not a number: the contact is face to face and the
resultant migrates as the faces slide apart. Retention at µ = 0.20 is 7.23 to 8.46 N. The
effective angles move less than a degree across that bracket, so the finding does not
depend on where the resultant is assumed to act.

**Runs, nine of them.** `lance_only` (prescribed tip lift — must return Stage 1 plus
1.06 %); `extract_mu000`, `extract_mu020`, `extract_mu030`; `insert_mu020` (full stroke
through the snap, on a 30 ms ramp because the stroke is four times longer and speed is
what the energy check sees); `extract_tpa_g010` (blocked); `extract_tpa_g085` (**not**
blocked — 0.85 mm is above the 0.757 mm limit derived below, and the same run doubles as
the control, since a contact that is defined but never reached has to cost exactly
nothing); `extract_mu020_fine` (element over the tooth halved) and `extract_mu020_slow`
(ramp doubled), which are the two sensitivity runs.

The TPA pair states the design criterion as a run rather than as a sentence. The criterion
is *not* clearance < protrusion — that was the first version of it, and the solver
contradicted it. See **The TPA criterion the model corrected** below.

Each deck is written into its own folder, `ls-dyna/stage3/stage3_<case>/`. LS-DYNA writes
`glstat`, `spcforc`, `nodout` and `d3plot` under fixed names into the directory it is
started in, so nine decks in one folder would overwrite each other's results. The
extractor takes the parent folder and finds each case by folder name — exact name first,
because `extract_mu020` is a substring of `extract_mu020_fine` and a substring match alone
would hand the baseline the sensitivity run's output.

**Mesh.** The tooth breakpoints land exactly on mesh stations, so both faces are planes
rather than staircases. The generator reads the angles back out of the finished mesh and
refuses to write a deck whose faces are not at the angles they were drawn at.

**A deck that was card-perfect and did nothing.** The propped-cantilever decks ran to
termination and wrote zeros — no reaction force, no internal energy, no error, no warning.
A part-ID collision: Stage 1 numbers its tip nodal rigid body part 2 because Stage 1 has no
second part, and that was carried into a deck where part 2 was a rigid part, so
`*BOUNDARY_PRESCRIBED_MOTION_RIGID` pointed at something `*MAT_RIGID` had locked in all six
directions. A nodal rigid body occupies a part ID like any other part. The audit had
compared keyword counts against the working decks and found them consistent, which they
were — what was wrong was what the IDs pointed at, and card counting cannot see that.
**And a second one of the same shape.** With the locking cycle built, `extract_mu020` also
ran to termination doing nothing — the terminal never moved, so it never touched the lance.
`*MAT_RIGID`'s CON1 field with CMO = 1.0 is a **code, not a bitmask**:

| CON1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| locks | none | x | y | z | x+y | **y+z** | z+x | x+y+z |

The terminal is driven along x and needs **5**. It was written as **6**, read as "everything
except x"; 6 locks z *and x*, so the prescribed motion had nothing left to move. Stage 2's
plate is driven the same way and had 5 in it the whole time — the value was on screen and
was not copied.

`check_ids()` now reads the finished deck text back and refuses a part-ID collision,
prescribed motion on the deformable part or on an undefined curve, a contact naming a part
that does not exist, **or a rigid part driven along a direction its own material card has
already constrained**. Each check was written after the failure it describes and tested by
reintroducing it.

**The common thread in both.** Neither failure produced an error, a warning or a zero-length
output file. Both produced a complete, clean run full of zeros. A solver that terminates
normally is not evidence that the model is doing anything, and neither is a deck whose cards
are all individually correct — what has to be checked is what the fields refer to.

**And a third failure, of a different kind — this one had a mechanism.** `insert_mu020`
terminated at 10.1 ms of a 20 ms ramp on excessive element distortion, with the tooth
crushed locally. `*CONTACT_AUTOMATIC_ONE_WAY_SURFACE_TO_SURFACE` checks slave **nodes**
against master **segments**, and nothing else. The terminal's leading edge is a sharp rigid
corner that travels the entire length of the tooth — and a master corner can sit inside a
slave element face, between its nodes, completely undetected. It had been gouging since
about 4 ms; the run only stopped once an element finally inverted, by which time the corner
was past the crest and under the retention face.

Unlike the first two, this one is a real modelling decision rather than a typo, and the fix
is four things:

| change | why |
|---|---|
| `*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE` — two-way | the terminal's corner is now checked as a node against the lance's faces |
| 0.15 mm × 45° chamfer on the terminal's leading edge | a real terminal has one; it turns a 90° corner into two 135° ones. Kept short and steep on purpose, so the tooth's 30° face stays the shallower surface and remains the angle the insertion force is predicted from |
| tooth mesh refined to 0.10 mm | an undetected corner penetrates by about one element length before a node notices |
| insertion ramp 20 → 30 ms | peak speed was 352 mm/s against Stage 2's 212; it is now 244 |

Only the first is strictly necessary. The others reduce how hard the contact has to work,
and the chamfer is there because the part has one.

**Why Stage 2 never hit this.** There the inclined face was on the rigid master and the
deformable tip rode it as a slave node on a large master face — the robust arrangement.
Stage 3 inverts it deliberately, because on a real lance the ramp is on the lance. Putting
the physics back where it belongs also put the contact into its fragile configuration, and
one-way contact was no longer good enough for it.

**A fourth failure, and this one is a modelling-limits result rather than a mistake in a
card.** The two TPA runs distorted: the lance was dragged along the insertion axis and
stretched. The cause is that **the terminal is rigid and kinematically driven, so once the
lance is resting on the TPA there is nothing in the model able to stop the terminal.** It
ploughs on. The blocked run had been driven 0.90 mm against a stop the lance reaches at
0.079 mm — eleven times past it.

That is the same objection that killed the very first Stage 3 plan, reappearing in a
different form. A prescribed displacement through a blocked path delivers unbounded force
just as a prescribed force through a blocked path delivers unbounded travel. A blocked run
is now driven to the stop plus 0.05 mm and no further. The overrun is not padding: it is
where the force rises steeply and the load path moves out of the lance's bending and into
the TPA, which is the thing the stage set out to measure. Past it the linear-elastic lance
has fractured and the model is describing a part that no longer exists.

#### The TPA criterion the model corrected

**And the FE model corrected the design criterion.** The control run was meant to be a TPA
fitted and useless — 0.65 mm of clearance against a 0.60 mm tooth protrusion, so the lance
should release with the TPA in place. It blocked. The one-line criterion *clearance <
protrusion* compares the clearance against how far the **tooth** must move, and ignores
that the TPA sits somewhere else on a beam that bends. Under the extraction load the lance
tip lifts **1.26 times** the crest, so a TPA reaching the tip stops the lance at 0.757 mm
of clearance, not 0.600.

| TPA reaches | blocks below |
|---|---|
| the crest | 0.600 mm |
| the lance tip | **0.757 mm** |

Corrected criterion: **the clearance must be compared against the lift at the most-lifted
station the TPA actually covers.** That is a statement about where the TPA is placed, not
about the tooth — and putting it over the tip buys 26 % more clearance for the same
function, which is a design lever the first version hid. The control run is now 0.85 mm,
and the generator refuses any TPA clearance within 10 % of the blocking limit, because a
run that close proves neither outcome.

**A fifth failure, and the only one that produced numbers rather than a crash.** The
first complete set of nine ran cleanly by every screen this document declares: energy
balance 0.99983 to 1.00000, hourglass exactly zero, all nine finished, the TPA control
matched the baseline to every decimal printed. Every force number in it was still void.

All output was written every 10 µs. The contact rings far above 50 kHz, so the force
history aliased. Aliased data cannot be filtered clean, because the folded content is
already inside the passband.

Two independent signatures identified it, and both belong in the screening list for
any explicit contact result:

1. **The answer moved with the filter.** The effective angle read 41.5° through a
   narrow moving average and 29.0° through a wide one, a 12° swing, while the same
   angle taken from the deflected shape moved 0.9° across the same range. A quantity
   that depends on the post-processing is not a measurement.
2. **The apparent frequency moved with the sample rate.** 44 kHz in one run, 16 kHz in
   the same model run at half the loading rate. A structural mode cannot do that. Only
   an alias can.

**Corrective action.** `*DATABASE_SPCFORC` and `*DATABASE_RCFORC` are written every
0.5 µs (2 MHz); `*DATABASE_NODOUT`, `*DATABASE_GLSTAT` and the rest stay at the previous
rate, because displacements are the double integral of the acceleration and the ring is
small in them. 1 µs was tried first and left some runs at 7 samples per cycle.

Resolving the ring was necessary and not sufficient. The ripple is 55 to 72 % of the
force's own level, so W/P read at an instant scatters by degrees however it is filtered.
The force is therefore read the way Stage 2 reads it — the ratio of the mean W to the
mean P over a 150 µs window — at five stations placed by the lance's lift, each compared
against the rotation over the same window. Forces are low-pass filtered at 10 kHz (SAE
J211 zero-phase Butterworth) before the window, to take the ring out.

**Three gates are applied to every run before a force may be quoted:**

| Gate | Threshold | Rationale |
|---|---|---|
| Samples per cycle of the ring | ≥ 8, where the ripple exceeds 5 % of the level | below this a low-pass has no usable transition band under Nyquist |
| Movement across averaging windows of 50–400 µs | ≤ 1.5° | more than this and the window is answering, not the model |
| Force route against shape route, worst station | ≤ 1.5° | the two share no data, so this is the error bar |

The first threshold was fixed before any run it was applied to. The other two were set
when the windowed estimator was introduced, after a first look at the 0.5 µs data, and
were not changed afterwards — including when `insert_mu020` came in just outside them.
A run that fails any gate has its force columns withheld. Its displacement measurements
are unaffected and are still reported.

**Two further errors in the same set.** The regression run `lance_only` was compared
against beam theory and read +3.8 %; the correct baseline is Stage 1's own FE result at
the same mesh — Stage 1 measured and recorded a 4.17 % model-form gap to Timoshenko, and
charging Stage 3 for it counts it twice. Against Stage 1 FE × 1.0106 for the tooth, the
run lands at −0.34 %. And `extract_mu020_fine` was read at the wrong node locations: one
node map was written for all nine cases, and that case builds its own mesh. Its forces
were unaffected because the root nodes keep their IDs, which is what made it silent.

**Verification loads, not service loads.** The Section 7 limitation applies to this stage
as to the others: the linear-elastic model bounds the forces from above and says nothing
about the load at which a real PBT-GF30 lance would fail.

---

## 6. Verification plan

| Check | Where | Criterion |
|---|---|---|
| Prescribed motion delivered | Stage 1 | tip node reaches −0.600 mm exactly |
| Closed-form agreement | Stage 1 | within ±2 % of P_T |
| Deflected shape | Stage 1 | quarter-point deflections match the closed form |
| Mesh convergence | Stage 1 | force stable with refinement |
| Element formulation | Stage 1 | ELFORM −1 on target; ELFORM 2 stiffens |
| Equilibrium | Stages 1–3 | `SPCFORC` = `BNDOUT` |
| Wedge relation | Stage 2 | W/P within ±2 % of (tan α + µ)/(1 − µ tan α), read while sliding |
| Contact formulation | Stage 2 | µ = 0: W/P = tan α, pure geometry |
| **Stage 1 regression** | Stage 2 | **µ = 0, β = 45° gives W/P = 1.000** |
| Sliding interface energy | Stage 2, µ = 0 | < 5 % of internal — penalty energy only |
| Frictional work | Stage 2, µ > 0 | sliding energy matches ½ µ N s / cos α |
| Energy balance | Stage 2 | total / initial energy within ±1 % |
| Contact vs constraint forces | Stage 2 | `RCFORC` = `SPCFORC` in x and z |
| Kinetic / internal energy | all explicit | small — quasi-static confirmed |
| Hourglass / internal energy | all explicit | < 10 % |
| Total energy | all explicit | flat — no contact energy injection |
| Contact-region mesh | Stage 2 | force not a mesh artefact |
| **Limiting case** | Stage 3 | `lance_only` reproduces Stage 1 plus the 1.1 % the tooth adds |
| Frictionless geometry | Stage 3 | W/P at µ = 0 returns tan of the *effective* angle |
| Retention and insertion | Stage 3 | inside the Stage 1 × Stage 2 bracket |
| Rotation correction | Stage 3 | measured face rotation matches the beam solution |
| TPA criterion | Stage 3 | blocks at 0.10 mm clearance, does not at 0.85 mm |

**Outcome of the closed-form criterion, recorded rather than revised.** The ±2 % band was
set against the Timoshenko value before any model was run, and the finest mesh sits at
+3.71 %, so **the criterion was not met**. The band stays on record. It is not used as a
pass/fail measure in the result figures, because it compares a 3D geometrically
non-linear solid against a 1D small-displacement solution, and the model carries a finite
width and a 6.4° tip rotation that the reference cannot represent. What is reported
instead, and separately: the convergence behaviour, the formulation sensitivity (0.75 %),
the discretization uncertainty (0.55 % at an assumed order), and the deviation from the
reference decomposed into named contributions. Numerical uncertainty and model form are
different quantities and are not summed into one number.

**Outcome of the Stage 3 checks, recorded the same way.**

| Check | Outcome |
|---|---|
| Limiting case | met — `lance_only` −0.34 % against Stage 1 FE × the tooth stiffening |
| Frictionless geometry | **not met** — `extract_mu000` fails the force gates; with µ = 0 there is no damping in the model and the contact rings undamped |
| Retention and insertion inside the bracket | **not assessed** — the force is read over a window at 0.90 of the lift rather than at release, and the insertion force fails the gates |
| Rotation correction | measured 7.67–7.78° at the retention face against 7.27° from the beam solution, 6 % high; 7.35° against 7.78° at the lead-in, 6 % low. Effective angles land 0.4–0.5° below prediction on both faces |
| TPA criterion | met — blocked at 0.10 mm, released at 0.85 mm |

---

## 7. Limitations

**The linear-elastic Stage 1 model is a solver-verification model, not a validated
prediction of the physical retention force.** The design strain of 1.125 % lies in a
region where the real PBT-GF30 stress–strain response deviates significantly from the
initial tangent modulus: the datasheet implies a secant modulus at break of 4,567 MPa
against a 9,800 MPa initial tangent — the curve softens by roughly 53 %. At the design
strain the linear model places the root at 80 % of break stress while at only 38 % of
break strain, which a real curve cannot do.

**Agreement between the analytical and FE Stage 1 forces therefore demonstrates
consistency of the numerical implementation under the linear-elastic assumption. It does
not establish the absolute force of a physical connector.** All forces quoted are upper
bounds.

Further limitations:

- The root face is clamped in all three translations, which also prevents the Poisson
  contraction of the root cross-section that beam theory allows. This was tested by
  releasing the lateral degrees of freedom on the root face and produced no measurable
  change in the tip force, so it is not a significant contributor here and the simple
  clamp is kept. Recorded because the test was run, not because the effect is large.
- The section has finite width. At b/t = 3.1 it sits between narrow-beam behaviour
  (modulus E) and plate behaviour (E/(1−ν²), +14 % at ν = 0.35), so some difference from
  a one-dimensional theory is expected and **does not reduce with mesh refinement**. It
  is a modelling difference, not a discretization error, and the two are not reported as
  one number.
- The comparison is against small-displacement beam theory, while the explicit solution
  is geometrically non-linear throughout. At δ/L = 0.075 the beam genuinely stiffens by
  of order 1 %. Real, expected, and not an error in either.
- Friction coefficient is assumed, not measured. Reported as a sensitivity band.
- Tensile modulus is used in a bending formulation. Flexural modulus is not published
  for this grade, and a thin moulded glass-filled section has skin–core fibre
  orientation that alters the effective bending stiffness.
- Long-term stiffness loss is indicated by the datasheet creep modulus but not modelled.
  Stress relaxation is the strictly correct mechanism and is not simulated.
- Geometry is simplified and representative, not a production part.
- **No experimental correlation.** The work is *verified*, not *validated*. Validation
  requires measurement — friction coefficient first, then a physical pull-out test.

---

## 8. Open items

1. Cite the specific design guide behind the 0.50 permissible-strain factor.
2. If buffer days permit: fit a two-point flow curve from the datasheet, run one case
   with `MAT_024`, and quantify how far the linear-elastic force is an overestimate.
   This converts the principal limitation into a result.

---

*Geometry synthetic and representative. Educational and self-learning use.*
