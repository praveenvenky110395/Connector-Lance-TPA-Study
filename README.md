# Connector Lance & TPA Study

A staged FE verification study of a plastic connector locking system — terminal, locking
lance, and secondary lock (TPA) — built so that **every added physical effect is verified
before the next one is introduced**.

**Status: Day 0.** Analytical reference and numerical specification only. No FE model
exists yet. That is deliberate — see *Why the analytical work comes first*, below.

---

## The engineering question

A connector terminal is held in its cavity by a thin flexible plastic finger: the
**lance**. The lance is weak, and the material's stiffness degrades over time. So a
second part — the **TPA** (Terminal Position Assurance) — is inserted afterwards to block
the lance from deflecting, and the terminal can no longer escape.

> **How does a connector locking feature behave when it stops being a simple bending
> problem and becomes a contact and friction problem — and what does the secondary lock
> actually change about the load path?**

---

## Method

Three stages. Each adds exactly one new piece of physics, and each is checked against
something independent before the next is added.

| Stage | Physics added | Checked against |
|---|---|---|
| **1 · Structural verification** | bending only, no contact | closed-form cantilever solution, plus a second solver |
| **2 · Contact mechanics** | contact, then friction | wedge relation; µ = 0 removes the friction term |
| **3 · Functional retention** | the TPA | no closed form — numerical verification and controlled comparison |

Stage 3 has no analytical solution. That is exactly why Stages 1 and 2 exist: by the time
the model is used where nothing can check it, it has already been checked twice.

### Why the analytical work comes first

The analytical reference is committed **before any FE model is built**, so the expected
values are a prediction rather than a description of results already obtained. Two
independent solvers agreeing would only prove they share an assumption; the closed form is
independent of both.

The prediction is specific: FE with a perfect root clamp should reproduce **3.6455 N**,
about **0.8 % below** the Euler–Bernoulli value of 3.6750 N, because of transverse shear.
Predicting the direction and size of the expected deviation is a different claim from
hoping the result lands close.

---

## Analytical targets

Computed by [`scripts/analytical.py`](scripts/analytical.py). Nothing is hard-coded; the script re-derives every
quantity by a second independent route and raises if the two disagree.

| Quantity | Value |
|---|---|
| Euler–Bernoulli tip force | 3.6750 N |
| **Timoshenko tip force — the FE target** | **3.6455 N** |
| Shear correction (κ = 5/6) | 0.8100 % |
| Root strain | 1.1250 % (75 % of permissible) |
| First bending mode | 5,110 Hz |
| Insertion force, µ = 0 | 2.1218 N |
| Retention force, µ = 0 | 3.6750 N — returns P exactly |

```
python scripts/analytical.py
```

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
tangent — the stress–strain curve softens by roughly 53 %. At the design strain of
1.125 % the linear model places the root at 80 % of break stress while at only 38 % of
break strain, which no real curve does. **All forces quoted here are upper bounds.**

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
│   └── analytical.py            analytical reference — self-verifying
└── results/
    └── analytical_targets.json  machine-readable targets for downstream comparison
```

Coming in Stages 1–3: `abaqus/`, `ls-dyna/`, `results/figures/`, `docs/technical_report.pdf`.

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
