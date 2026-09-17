# Connector Lance & TPA Study

A staged FE verification study of a plastic connector locking system — terminal, locking
lance, and secondary lock (TPA) — built so that **every added physical effect is verified
before the next one is introduced**.

**Status: Stage 1 complete.** Analytical reference frozen, boundary conditions verified,
six runs done, results below and in [`results/`](results/).

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
| **2 · Contact mechanics** | contact, then friction | wedge relation; µ = 0 removes the friction term |
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

Computed by [`scripts/analytical.py`](scripts/analytical.py). Nothing is hard-coded; the script re-derives every
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
│   ├── make_stage1_models.py    writes the six LS-DYNA decks
│   └── postprocess_stage1.py    convergence order, Richardson, GCI, figures
├── ls-dyna/stage1/              six decks: four mesh levels, two extra formulations
└── results/
    ├── analytical_targets.json  machine-readable targets
    ├── stage1_results.csv       filled in from the runs
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
