# Numerical Specification

**Connector Lance & TPA Study** — frozen inputs and analytical expectations.

| | |
|---|---|
| Version | 1.0 |
| Status | Inputs frozen; Stage 1 decks generated, evaluation pending |
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

Terminal added (2.50 × 2.00 × 12.00 mm, CuSn6-type, E = 115,000 MPa, ν = 0.34,
ρ = 8.90e-9 t/mm³). Housing added — the resulting additional softening is quantified as
a result, not treated as error.

- `*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE`, `SOFT=1` baseline, `SOFT=2` comparison
- `VDC = 10–20` to damp contact chatter
- **µ = 0 first**, then 0.10 / 0.20 / 0.30
- Insertion stroke 4 mm at 400 mm/s over 10 ms (≈ 51 × T₁)
- Extract `RCFORC` at the interface and `BNDOUT` at the driven node

### Stage 3 — TPA

TPA block 2.50 × 0.80 × 4.00 mm, seated with **0.10 mm clearance** above the undeflected
lance. The lance can therefore lift 0.10 mm against the 0.60 mm it needs to release —
**83 % of the required deflection is blocked.**

- Case A: TPA off. Pull-out 1.5 mm at 300 mm/s.
- Case B: TPA on. Identical in every other respect.
- **Limiting-case regression:** with the TPA contact deactivated, the Stage 3 model must
  reproduce Stage 2. This is the primary verification argument for Stage 3.
- If time permits: correctly seated vs under-seated terminal, and TPA closure blocked by
  an under-seated terminal — the position-assurance function itself.

This is an **investigation, not a demonstration**. The question is how the TPA changes
the load path, not whether it improves retention. A smaller-than-expected effect, or a
shift in failure mode, is a valid and more interesting result.

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
| Wedge relation | Stage 2 | µ = 0 matches P·tan α |
| Stage 1 regression | Stage 2 | µ = 0, β = 45° returns P_EB |
| Kinetic / internal energy | all explicit | small — quasi-static confirmed |
| Hourglass / internal energy | all explicit | < 10 % |
| Total energy | all explicit | flat — no contact energy injection |
| Contact-region mesh | Stage 2 | force not a mesh artefact |
| **Limiting case** | Stage 3 | TPA deactivated reproduces Stage 2 |

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
