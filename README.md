# Connector Lance & TPA Study

A staged FE verification study of a plastic connector locking system — terminal, locking
lance, and secondary lock (TPA) — built so that **every added physical effect is verified
before the next one is introduced**.

**Status: Stages 1 and 2 complete.** Stage 2 matches the closed-form wedge relation to within
0.41 % across both ramp angles and a friction sweep to µ = 0.30. Analytical reference frozen,
boundary conditions verified, results below and in `results/`.

---

## The engineering question

A connector terminal is held in its cavity by a thin flexible plastic finger: the
**lance**. The lance is a flexible feature, and long-term material behaviour can reduce
its effective load-carrying capability. So a second part — the **TPA** (Terminal Position
Assurance) — is inserted afterwards to block the lance from deflecting, and the terminal
can no longer escape.

> **How does a connector locking feature behave when it stops being a simple bending
> problem and becomes a contact and friction problem — and what does the secondary lock
> actually change about the load path?**

---

## Method

Three stages. Each adds exactly one new piece of physics, and each is checked against
something independent before the next is added.

| Stage | Physics added | Checked against |
|---|---|---|
| **1 · Structural verification** | bending only, no contact | closed-form cantilever solution |
| **2 · Contact mechanics** | contact and friction | wedge relation, checked on W/P — a ratio with no stiffness in it |
| **3 · Functional retention** | the TPA | no closed form — numerical verification and controlled comparison |

Stage 3 has no analytical solution. That is exactly why Stages 1 and 2 exist: by the time
the model is used where nothing can check it, it has already been checked twice.

The analytical reference is committed **before any FE model is built**, so the expected
values are a prediction rather than a description of results already obtained. The closed
form is not a solver and cannot share a solver's mistakes. A second FE code was considered
and deliberately left out: it would have cost days of unfamiliar tooling to give a second
opinion on a question the analytical solution already answers.

---

## Analytical targets

Computed by `scripts/analytical.py`. Nothing is hard-coded; the script re-derives every
quantity by a second independent route and raises if the two disagree.

| Quantity | Value |
|---|---|
| Euler–Bernoulli tip force | 3.6750 N |
| **Timoshenko tip force — analytical reference** | **3.6455 N** |
| Shear correction (κ = 5/6) | 0.8100 % |
| Root strain | 1.1250 % (75 % of permissible) |
| First bending mode | 5,110 Hz |
| Insertion force, µ = 0 | 2.1218 N |
| Retention force, µ = 0 | 3.6750 N — returns P exactly |

```
python scripts/analytical.py
python scripts/make_stage1_models.py
python scripts/postprocess_stage1.py
```

---

## Stage 1 finding: how the tip deflection is applied

The obvious way to impose a tip deflection on a solid model is to prescribe the same z
displacement at every node of the tip face. On this geometry that is wrong, and the reason
is geometric rather than numerical.

At a tip deflection of 0.6 mm on an 8 mm cantilever, the tip section rotates **6.4°**. A
plane section rotated by θ projects shorter onto the z axis by t(1 − cos θ). Holding every
node of that face at the same z forbids the projection change, so the material has to
stretch instead:

> **ε_zz = θ² / 2 = 0.63 %**

against a root bending strain of 1.125 %. An artificial through-thickness strain worth
more than half the peak design strain of the model, injected at the tip, and paid for in
reaction force. On the n4 mesh it costs about **7 %**.

The model therefore ties the tip face into a rigid section and prescribes the deflection
at a control node on the neutral axis, with all rotations free. This is Euler–Bernoulli's
"plane sections remain plane" written as a constraint — the FE model is made to adopt the
assumption the closed form is built on, which is the point of a verification study.

The effect scales as θ², so it is invisible at small deflections and grows quickly. Here
δ/L = 0.075, which is where it stops being negligible.

**Two other explanations were tested first and set aside:**

| Hypothesis | Test | What the data shows |
|---|---|---|
| Element formulation — shear locking was the initial suspicion | ELFORM −1, 1 and 2 on a common mesh | 0.75 % spread. Formulation sensitivity is small, so the formulation choice is not driving the result. |
| Over-constrained root clamp | lateral DOF released on the root face, one node held to prevent drift | no measurable change in tip force |

---

## Stage 1 results

> **The mesh sequence settles at 3.7809 N, and the result is insensitive to element
> formulation within 0.75 %. The converged 3D solution sits +3.71 % above the idealized
> Timoshenko reference.**

Tip force at a prescribed deflection of 0.600 mm. The control node reads −0.600 mm on
every run, checked before any force was recorded.

| Mesh | h [mm] | Force [N] | Δ vs Timoshenko |
|---|---|---|---|
| n1 | 0.800 | 3.0285 | −16.93 % |
| n2 | 0.400 | 3.8086 | +4.47 % |
| n4 | 0.200 | 3.7975 | +4.17 % |
| **n8** | **0.100** | **3.7809** | **+3.71 %** |

n1 is outside the asymptotic range — one linear element through the thickness carries
constant transverse shear and cannot represent the parabolic distribution — and is
excluded from the fit.

Element formulation on the common mesh (n4):

| Formulation | Force [N] | Δ vs Timoshenko |
|---|---|---|
| ELFORM 1 — one-point, hourglass IHQ 6 / QM 1.0 | 3.7917 | +4.01 % |
| ELFORM −1 — fully integrated, poor-aspect-ratio form | 3.7975 | +4.17 % |
| ELFORM 2 — fully integrated | 3.8202 | +4.79 % |

**Maximum formulation spread: 0.75 %.** Under-integrated softest, fully integrated
stiffest — the ordering element theory predicts, and small enough that the formulation
choice is not driving the result.

**Solution quality.** Kinetic over internal energy peaks at 0.007 % during the ramp.
Hourglass energy on ELFORM 1 is 0.069 % of internal, against a 10 % screening criterion.
Reaction force drifts by at most 0.03 % across the hold, so the response has settled.

**Convergence.** The three-mesh Richardson fit does not hold: the differences grow under
refinement (−0.0111 N then −0.0166 N) instead of shrinking, so at least one mesh is not
demonstrably in the asymptotic range. No observed order and no extrapolation are reported.

> Estimated discretization uncertainty: **≈ 0.55 %** (GCI on the two finest meshes,
> **assuming p = 1**. The observed convergence order was not used, because the meshes were
> not demonstrably in the asymptotic range.)

A likely reason the sequence is not a clean power law: the tip constraint is applied on a
face whose own discretization changes with the mesh, so refinement is not varying
discretization alone.

### Numerical sensitivity and model form are separate quantities

| | |
|---|---|
| Element formulation spread | 0.75 % |
| Estimated discretization uncertainty (assumed p = 1) | ≈ 0.55 % |
| **Deviation from the analytical reference** | **3.71 %** |

The numerical sensitivity of the model is several times smaller than its distance from the
closed form. **What separates the two is model form, not numerics.**

The remaining deviation is consistent with geometric non-linearity, finite-width
three-dimensional effects, and residual discretization uncertainty. Approximate sensitivity
estimates indicate contributions of the order of **1 %**, **2 %** and **< 1 %**
respectively. These are interpretive estimates of the individual mechanisms, not an exact
error budget — the effects are not independent and are not claimed to sum to the total.

Two of them are checked against evidence outside the force comparison. The finite-width
term places the result **20.6 %** of the way from narrow-beam behaviour (modulus E) to
plate behaviour (E/(1−ν²)), which for b/t = 3.1 is the expected region. The geometric term
rests on a tip run-in of 0.027 mm predicted from the slope field, which matches the tip
x-displacement in `nodout` once the +0.045 mm rotation term is removed.

The closed form is therefore used as a bracket as well as a reference: a section of this
aspect ratio should lie between the narrow-beam value (3.6455 N) and the plate-modulus
value (4.1544 N). The finest mesh lies at 27 % of that interval.

**On the initial benchmark.** A ±2 % acceptance band (3.5726 – 3.7184 N) was defined
against the Timoshenko value before any model was run. **It was not met.** It is kept on
record rather than redrawn, and it is not used as a pass/fail criterion in the figures: it
tests a 3D geometrically non-linear solid against a 1D small-displacement solution, and
once the two are known not to be the same problem, the informative measures are the
convergence behaviour, the formulation sensitivity, and the deviation interpreted above.

### Conclusion

> Stage 1 demonstrates a numerically stable LS-DYNA solid model for the prescribed
> cantilever deformation. The solution converges toward approximately 3.78 N, with less
> than 1 % sensitivity to element formulation and low dynamic, hourglass and settling
> contributions. The converged 3D result differs by approximately 3.7 % from the idealized
> Timoshenko reference. This deviation is interpreted as a model-form difference associated
> with the three-dimensional geometry and non-linear deformation, rather than as numerical
> instability.

---

## Stage 2: contact, and why the check is a ratio

The Stage 1 lance, unchanged — same geometry, material, mesh and root clamp, generated by
importing the Stage 1 generator so it cannot drift. What is added is a rigid plate above the
tip with an inclined lower face, driven along −x. Advancing by s lowers the face at the tip
by s·tan α, so the plate deflects the lance with nothing prescribed on the lance at all.

The inclination sits on the **rigid** body on purpose. Contact is on the lance's tip edge,
so the contact normal is the plate's normal — and the plate does not rotate. Put the wedge
face on the lance and the 6.4° tip rotation from Stage 1 would tilt the normal with it; at
30° that is a 28 % error in tan α. The same rotation that broke Stage 1's boundary condition
would have broken Stage 2's reference.

The comparison is on a ratio:

> **W/P = (tan α + µ) / (1 − µ tan α)**

W is the summed axial reaction at the root, P the summed transverse one, both from `spcforc`
at the same instant. Global equilibrium makes them the two components of the contact force,
so the relation follows from the contact normal and the friction law alone.

**The ratio holds no stiffness** — not E, not I, not the deflection. The Stage 1 difference
from beam theory is a stiffness difference, so it cancels out of this comparison exactly.
Stage 1 can sit 3.7 % above beam theory and Stage 2 can still be exact; both are true at
once. That is what the staged structure was built for.

| α | µ | W/P target |
|---|---|---|
| 30° | 0 | 0.5774 |
| 30° | 0.10 / 0.20 / 0.30 | 0.7189 / 0.8788 / 1.0611 |
| **45°** | **0** | **1.0000** |
| 45° | 0.20 | 1.5000 |

The 45° frictionless case is the regression test. Nothing but the contact normal sets it —
no friction, no stiffness, no material. It gets run first: if it is wrong, the contact is
wrong and no other number in the stage is safe.

### Stage 2 results

| case | α | µ | W [N] | P [N] | W/P | closed form | deviation |
|---|---|---|---|---|---|---|---|
| a30_mu000 | 30° | 0 | 1.4847 | 2.5753 | 0.5765 | 0.5774 | −0.14 % |
| a30_mu010 | 30° | 0.10 | 1.8653 | 2.5843 | 0.7218 | 0.7189 | +0.41 % |
| a30_mu020 | 30° | 0.20 | 2.2648 | 2.5777 | 0.8786 | 0.8788 | −0.02 % |
| a30_mu030 | 30° | 0.30 | 2.7411 | 2.5836 | 1.0610 | 1.0611 | −0.02 % |
| **a45_mu000** | **45°** | **0** | **2.6253** | **2.6249** | **1.0002** | **1.0000** | **+0.02 %** |
| a45_mu020 | 45° | 0.20 | 3.9279 | 2.6246 | 1.4966 | 1.5000 | −0.23 % |

Worst deviation across the six, **+0.41 %**, against a ±2 % band. KE/IE at most 1.01 %,
energy ratio 1.00000 on every run.

Getting there required finding out that the first measurement window was the wrong one, and
the way it was found is more useful than the numbers.

Read during the hold — the 0.5 ms after the plate stops, where the force is flat and the
obvious place to measure — the two frictionless runs came out at +0.03 % and −0.07 %, and
the four friction runs at **−0.89, −1.52, −1.90 and −0.57 %**. All inside the band, so by the
written criterion the stage passed. But the miss grew monotonically with µ and sat at roughly
four standard errors on each of the four, while the frictionless pair sat at a tenth of one.
A deviation that tracks the parameter being swept is not scatter.

The energy check then pointed the other way. Frictional work over the stroke came out 7–12 %
**above** the closed form, not below. Under-mobilised friction would have been short on both.
So friction was fully developed while the plate was moving, and the shortfall was somewhere
else.

W/P = (tan α + µ)/(1 − µ tan α) is a *sliding* relation — it comes from Coulomb friction with
relative motion at the interface. During the hold the plate is stationary, so the contact is
stuck, and a stuck penalty interface carries whatever tangential force its stick spring held
when motion stopped, less whatever it sheds as the lance settles back. With µ = 0 there is no
tangential force to relax, which is why exactly those two runs were unaffected.

The test was stated before the data was looked at: read while the plate is moving, the
friction runs should sit on the closed form and the µ-trend should disappear; if they stayed
2 % low in every window, the friction model itself was wrong and the hold number would stand.
They sit on the closed form. Over all 35 windows between 0.8 and 1.5 ms long inside 2.2–4.2 ms,
every case has a median deviation within 0.2 % and no single window anywhere exceeds 0.93 %,
so the answer does not depend on where in the stroke it is read — only on whether the plate
is moving.

![W/P through the stroke](results/figures/stage2_ratio_history.png)

Two other things in the data agree. Tip deflection falls as µ rises — 0.6014, 0.6008, 0.6000,
0.5990 mm — because friction drags the tip back along −x so it rides higher on the inclined
face. And the ratio scatter collapses when friction is switched on, 7.3 % to 1.1 % at 45°,
because friction damps the axial mode that was making the frictionless ratio noisy.

Extraction is scripted rather than clicked: `scripts/extract_stage2.py` reads `spcforc`,
`nodout` and `glstat` out of the six run folders and writes `results/stage2_results.csv` plus
the full reaction history. `spcforc` already closes each output block with a `force
resultants` line — the x, y, z sum over the 35 root nodes — so W and P are read from the file
rather than summed by hand. The script was written against the 45° frictionless run and
reproduces the numbers taken manually from LS-PrePost, 3.7302 N and 3.7329 N, to four
figures. Having the history in a file rather than in a GUI is what made the window question
answerable at all.

Reading the first run also corrected two acceptance criteria, both written for a Stage 1
without contact. The 5 % screen on sliding interface energy only means something with µ = 0,
where the counter holds penalty energy alone (1.90 % here); with friction on it also collects
real frictional work, which at 30° and µ = 0.3 is near 80 % of the internal energy. Those
cases are compared against the closed-form friction work instead, and the screen that holds
for every case is the energy balance, ±1 %. KE/IE is likewise read over the second half of
the ramp only — earlier the lance holds almost no internal energy and the ratio is set by its
denominator, reading 7.3 % at 1.02 ms and 0.32 % at 2.51 ms on the same run while the kinetic
energy is still rising.

All three changes are to how the result is measured and judged. Nothing in the model was
altered to meet a criterion, and both the hold and the sliding numbers are reported.

`VDC = 0` throughout — no contact damping. Damping would add force to the quantity being
measured, and if a run needs it to stay stable, that is a result to report rather than a
setting to bury.

---

## Geometry and material

Lance L = 8.00 × b = 2.50 × t = 0.80 mm, deflection y = 0.60 mm, L/t = 10.
Lead angle 30°, return angle 45°, **both measured from the insertion axis**.

Material: **BASF Ultradur® B 4300 G6 (PBT-GF30)** — E = 9,800 MPa (ISO 527),
strain at break 3.0 % (ISO 527), density 1,530 kg/m³ (ISO 1183).

Full specification, including the angle convention, the geometric consistency constraint
and the full limitations section: **[`docs/numerical_spec.md`](docs/numerical_spec.md)**.

---

## Principal limitation — stated up front

**The linear-elastic Stage 1 model is a solver-verification model, not a validated
prediction of a physical retention force.**

The datasheet implies a secant modulus at break of 4,567 MPa against a 9,800 MPa initial
tangent — the stress–strain curve softens by roughly 53 %. At the design strain of 1.125 %
the linear model places the root at 80 % of break stress while at only 38 % of break
strain. The forces quoted are analytical reference values under the stated linear-elastic
assumptions.

Agreement between the analytical and FE forces demonstrates consistency of the numerical
implementation under the linear-elastic assumption. It does not establish the absolute
force of a physical connector.

This work is **verified, not validated**. Validation requires measurement.

---

## Repository layout

```
├── README.md
├── docs/
│   └── numerical_spec.md        frozen inputs, targets, verification plan, limitations
├── scripts/
│   ├── analytical.py            analytical reference — self-verifying
│   ├── make_stage1_models.py    writes the six Stage 1 decks
│   ├── postprocess_stage1.py    convergence order, GCI, figures
│   ├── make_stage2_models.py    writes the six Stage 2 decks
│   ├── extract_stage2.py        reads the run folders, fills stage2_results.csv
│   └── postprocess_stage2.py    W/P against the wedge relation, friction figure
├── ls-dyna/
│   ├── stage1/                  four mesh levels, two extra formulations
│   └── stage2/                  two ramp angles, friction sweep
└── results/
    ├── analytical_targets.json  machine-readable targets
    ├── stage1_results.csv       filled in from the runs
    ├── stage2_results.csv
    └── figures/
```

**Tools: LS-DYNA (Student) and Python.** Deliberately kept to two, both of which every
result here can be fully defended in.

---

## Related work

Preceded by **[LS-DYNA Crimping Simulation](https://github.com/praveenvenky110395/LS-DYNA-Crimping-Simulation)**
— crimp forming of a 19-strand conductor. That study covered *manufacturing*: how the
terminal is made. This one covers *function*: how the terminal is held.

---

## Scope

Geometry is synthetic and representative of automotive connector locking features. It is
**not** a reproduction of any manufacturer's product, and no manufacturer data beyond
published material datasheets has been used. Educational and self-learning use.
