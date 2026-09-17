#!/usr/bin/env python3
"""
analytical.py - closed-form reference for the connector lance study.

Computes the analytical reference values for the locking lance ahead of the FE
model build. All values derive from the inputs block below; none are hard-coded.

Each quantity is computed by its compact closed form and re-derived by an
independent algebraic route; disagreement raises. The wedge relation is checked
against its two limiting cases.

Units: mm - tonne - s - N - MPa
    length      mm
    mass        tonne (Mg)
    time        s
    force       N
    stress      MPa = N/mm^2
    density     tonne/mm^3   (1 g/cm^3 = 1e-9 t/mm^3)
    frequency   Hz

Scope
    Stage 1 is linear-elastic with idealized boundary conditions. The absolute
    forces reported here are analytical reference values under the stated
    linear-elastic assumptions. No experimental validation exists; the project
    provides analytical and numerical verification only. See
    docs/numerical_spec.md, section "Limitations".

Author: Praveen Venkatesh Sethumadhavan Vasan
Use: educational and self-learning. Geometry is synthetic and representative; it
is not a reproduction of any manufacturer's product.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

# =============================================================================
# Inputs. All downstream values derive from this block.
# =============================================================================


@dataclass(frozen=True)
class Geometry:
    """Locking lance, terminal notch and ramp angles. All lengths in mm."""

    L: float = 8.00      # lance free length (root clamp to tooth centreline)
    b: float = 2.50      # lance width
    t: float = 0.80      # lance thickness
    y: float = 0.60      # prescribed deflection = tooth protrusion = notch depth

    # Ramp angles, both measured from the insertion axis. A shallow ramp is a
    # small angle: alpha -> 0 is a face parallel to the insertion direction
    # (sliding only), alpha -> 90 deg is normal to it (self-locking). The
    # opposite convention inverts the Stage 2 reference values.
    alpha_deg: float = 30.0   # lead-in ramp, insertion side
    beta_deg: float = 45.0    # return face, retention side

    @property
    def slenderness(self) -> float:
        """L/t. Screening criterion for beam-theory applicability."""
        return self.L / self.t


@dataclass(frozen=True)
class Material:
    """BASF Ultradur B 4300 G6, PBT-GF30. Datasheet values with ISO methods."""

    name: str = "BASF Ultradur B 4300 G6 (PBT-GF30)"
    E: float = 9800.0          # MPa, tensile modulus, ISO 527-1/-2
    nu: float = 0.35           # assumed; not on the datasheet
    rho: float = 1.53e-9       # tonne/mm^3, from 1530 kg/m^3, ISO 1183
    sigma_break: float = 137.0  # MPa, stress at break, ISO 527-1/-2
    eps_break: float = 0.030   # strain at break, ISO 527-1/-2
    E_creep_1000h: float = 7500.0  # MPa, tensile creep modulus 1000 h/23C, ISO 899-1

    # Permissible strain for a one-time-assembly snap-fit in a glass-reinforced
    # grade, as a fraction of the datasheet strain at break. 0.50 is the
    # conservative end of the published design-guide range.
    # OPEN ITEM: cite the specific design guide consulted.
    permissible_strain_factor: float = 0.50

    @property
    def G(self) -> float:
        """Shear modulus from isotropic elasticity, MPa."""
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def eps_permissible(self) -> float:
        return self.permissible_strain_factor * self.eps_break

    @property
    def E_secant_at_break(self) -> float:
        """Secant modulus at the datasheet break point, MPa."""
        return self.sigma_break / self.eps_break


# Shear correction factor for a rectangular cross-section.
#   "classic" : kappa = 5/6
#   "cowper"  : kappa = 10(1+nu)/(12+11nu)
# Difference in the resulting reference force here is ~0.02 %.
KAPPA_CONVENTION = "classic"


# =============================================================================
# Closed-form relations
# =============================================================================


def second_moment_of_area(g: Geometry) -> float:
    """I = b t^3 / 12  [mm^4]  - rectangular section about the bending axis."""
    return g.b * g.t**3 / 12.0


def shear_coefficient(m: Material, convention: str = KAPPA_CONVENTION) -> float:
    """Timoshenko shear correction factor kappa, rectangular section."""
    if convention == "classic":
        return 5.0 / 6.0
    if convention == "cowper":
        return 10.0 * (1.0 + m.nu) / (12.0 + 11.0 * m.nu)
    raise ValueError(f"unknown kappa convention: {convention!r}")


def tip_force_euler_bernoulli(g: Geometry, m: Material) -> float:
    """
    Tip force for prescribed deflection y.  [N]

        delta = P L^3 / (3 E I)   ->   P = E b t^3 y / (4 L^3)

    Bending only: no shear deformation, rigid root, tip point load.
    """
    return m.E * g.b * g.t**3 * g.y / (4.0 * g.L**3)


def shear_correction_ratio(g: Geometry, m: Material) -> float:
    """
    Shear contribution to tip deflection, as a fraction of the bending
    contribution:

        delta_shear / delta_bending = 3 E I / (kappa G A L^2)

    For a rectangular section with kappa = 5/6 this reduces to

        = 0.6 (1 + nu) (t/L)^2
    """
    I = second_moment_of_area(g)
    A = g.b * g.t
    kappa = shear_coefficient(m)
    return 3.0 * m.E * I / (kappa * m.G * A * g.L**2)


def tip_force_timoshenko(g: Geometry, m: Material) -> float:
    """
    Tip force including transverse shear, same prescribed deflection.  [N]

    Shear deformation lowers the required force relative to Euler-Bernoulli.
    This is the analytical reference for a well-resolved FE model with the same
    idealized boundary conditions.
    """
    return tip_force_euler_bernoulli(g, m) / (1.0 + shear_correction_ratio(g, m))


def deflected_shape(xi: float, shear_ratio: float) -> float:
    """
    Deflection at x, divided by the deflection at the tip.  [-]

        w(x)/w(L) = [ (3 xi^2 - xi^3) + 2 r xi ] / [ 2 (1 + r) ],    xi = x/L

    The cubic term is bending. The term in r is transverse shear, which varies
    linearly along the beam because the shear force is constant under a tip load.

    Dividing by the tip value makes the function return 1 at xi = 1 whatever r
    is, so the shape can be compared against the FE model without the force
    entering the comparison at all. The tip deflection is prescribed in the FE
    model, so this is a check on the response between the root and the tip -
    independent of the reaction force, and therefore not a restatement of it.

    Pass shear_ratio = 0 for the Euler-Bernoulli shape.
    """
    return ((3.0 * xi**2 - xi**3) + 2.0 * shear_ratio * xi) / (2.0 * (1.0 + shear_ratio))


def root_stress(g: Geometry, m: Material) -> float:
    """Peak bending stress at the root, sigma = M c / I.  [MPa]"""
    P = tip_force_euler_bernoulli(g, m)
    return (P * g.L) * (g.t / 2.0) / second_moment_of_area(g)


def max_root_strain(g: Geometry) -> float:
    """
    Peak bending strain at the root.  [-]

        eps = 3 t y / (2 L^2)

    Independent of E; compared against the permissible strain.
    """
    return 3.0 * g.t * g.y / (2.0 * g.L**2)


def wave_speed(m: Material) -> float:
    """One-dimensional elastic wave speed, c = sqrt(E/rho).  [mm/s]"""
    return math.sqrt(m.E / m.rho)


def first_bending_frequency(g: Geometry, m: Material) -> float:
    """
    First cantilever bending mode.  [Hz]

        f1 = (lambda1^2 / 2 pi) sqrt( E I / (rho A L^4) ),   lambda1 = 1.8751

    Explicit time integration is used for the FE model; quasi-static response is
    assessed from the loading duration relative to T1 = 1/f1 and from the energy
    balance.
    """
    lam1 = 1.8751040687119611
    I = second_moment_of_area(g)
    A = g.b * g.t
    return (lam1**2 / (2.0 * math.pi)) * math.sqrt(m.E * I / (m.rho * A * g.L**4))


def wedge_force(P: float, angle_deg: float, mu: float) -> float:
    """
    Axial force to drive a spring-loaded tooth along an inclined face.  [N]

        W = P (tan(angle) + mu) / (1 - mu tan(angle))

    Free-body derivation:
        transverse equilibrium :  N (cos a - mu sin a) = P
        axial force            :  W = N (sin a + mu cos a)

    Angle measured from the insertion axis:
        angle -> 0   : face parallel to insertion, W -> mu P
        angle -> 90  : face normal to insertion, W -> infinity (self-locking)

    Rigid-body relation with a single contact point and Coulomb friction. Used as
    an idealized analytical reference for Stage 2, not as a prediction of the 3D
    contact model.

    P is the transverse spring force from the lance at the deflection reached on
    the face.
    """
    ta = math.tan(math.radians(angle_deg))
    denom = 1.0 - mu * ta
    if denom <= 0.0:
        return math.inf  # self-locking
    return P * (ta + mu) / denom


# =============================================================================
# Self-checks - each result re-derived by an independent route
# =============================================================================


def _self_check(g: Geometry, m: Material, rtol: float = 1e-9) -> list[str]:
    """Re-derive each quantity a second way. Raises on disagreement."""
    checks: list[str] = []

    def close(a: float, b: float) -> bool:
        return abs(a - b) <= rtol * max(1.0, abs(b))

    # Tip force: compact form vs 3 E I y / L^3
    I = second_moment_of_area(g)
    p_compact = tip_force_euler_bernoulli(g, m)
    p_from_EI = 3.0 * m.E * I * g.y / g.L**3
    assert close(p_compact, p_from_EI), "tip force: forms disagree"
    checks.append("P_EB: E b t^3 y/(4L^3) == 3 E I y/L^3")

    # Strain: compact form vs sigma/E through M c / I
    eps_compact = max_root_strain(g)
    eps_from_stress = root_stress(g, m) / m.E
    assert close(eps_compact, eps_from_stress), "strain: forms disagree"
    checks.append("eps: 3ty/(2L^2) == (M c / I)/E")

    # Shear correction: general form vs the rectangular reduction
    if KAPPA_CONVENTION == "classic":
        general = shear_correction_ratio(g, m)
        simplified = 0.6 * (1.0 + m.nu) * (g.t / g.L) ** 2
        assert close(general, simplified), "shear correction: forms disagree"
        checks.append("shear: 3EI/(kappa G A L^2) == 0.6(1+nu)(t/L)^2")

    # Deflected shape: zero at the root and unity at the tip, for either theory
    r = shear_correction_ratio(g, m)
    assert close(deflected_shape(0.0, r), 0.0), "shape: not zero at the root"
    assert close(deflected_shape(1.0, r), 1.0), "shape: not unity at the tip"
    checks.append("shape: w(0)/w(L) = 0 and w(L)/w(L) = 1")

    # Shape must rise monotonically from root to tip
    shape = [deflected_shape(i / 20.0, r) for i in range(21)]
    assert all(a < z for a, z in zip(shape, shape[1:])), "shape: not monotonic"
    checks.append("shape: monotonic from root to tip")

    # Mid-span with no shear: the cubic gives exactly 5/16
    assert close(deflected_shape(0.5, 0.0), 5.0 / 16.0), "shape: mid-span value wrong"
    checks.append("shape: Euler-Bernoulli mid-span == 5/16 exactly")

    # Wedge limit: face parallel to insertion reduces to sliding friction
    P = p_compact
    assert close(wedge_force(P, 0.0, 0.25), 0.25 * P), "wedge: alpha=0 limit fails"
    checks.append("wedge limit: alpha=0 -> W = mu P")

    # Wedge limit: 45 deg face, frictionless, returns the transverse spring force
    assert close(wedge_force(P, 45.0, 0.0), P), "wedge: 45 deg frictionless fails"
    checks.append("wedge limit: beta=45 deg, mu=0 -> W = P  [Stage 2 regression test]")

    # Wedge force must increase monotonically with friction
    ws = [wedge_force(P, g.alpha_deg, mu) for mu in (0.0, 0.1, 0.2, 0.3)]
    assert all(x < z for x, z in zip(ws, ws[1:])), "wedge: not monotonic in mu"
    checks.append("wedge: monotonic in mu")

    # Slenderness screening criterion
    assert g.slenderness >= 8.0, \
        f"L/t = {g.slenderness:.1f} is below the slenderness screening criterion"
    checks.append(f"slenderness: L/t = {g.slenderness:.1f} >= 8")

    return checks


# =============================================================================
# Report
# =============================================================================


def build_targets(g: Geometry, m: Material) -> dict:
    """Analytical reference values in one dictionary, for downstream comparison."""
    P_EB = tip_force_euler_bernoulli(g, m)
    P_T = tip_force_timoshenko(g, m)
    f1 = first_bending_frequency(g, m)
    eps = max_root_strain(g)
    r = shear_correction_ratio(g, m)

    # Quarter points along the length. These are the stations the FE deflected
    # shape is read at; they exist as nodes on every mesh except n1, where the
    # element spacing is 0.8 mm and 2 mm and 6 mm fall between nodes.
    stations = (0.25, 0.50, 0.75, 1.00)

    return {
        "units": "mm-tonne-s-N-MPa",
        "kappa_convention": KAPPA_CONVENTION,
        "geometry": asdict(g),
        "material": asdict(m),
        "stage1": {
            "I_mm4": second_moment_of_area(g),
            "P_euler_bernoulli_N": P_EB,
            "P_timoshenko_N": P_T,
            "shear_correction_pct": r * 100.0,
            "root_stress_MPa": root_stress(g, m),
            "root_strain_pct": eps * 100.0,
            "permissible_strain_pct": m.eps_permissible * 100.0,
            "strain_utilisation_pct": eps / m.eps_permissible * 100.0,
            "f1_Hz": f1,
            "T1_ms": 1.0 / f1 * 1e3,
            "wave_speed_mm_s": wave_speed(m),
            "acceptance_band_N": [P_T * 0.98, P_T * 1.02],
            "deflected_shape_mm": {
                f"x_{xi * g.L:.1f}mm": {
                    "x_over_L": xi,
                    "w_euler_bernoulli": deflected_shape(xi, 0.0) * g.y,
                    "w_timoshenko": deflected_shape(xi, r) * g.y,
                }
                for xi in stations
            },
        },
        "stage2": {
            f"mu_{mu:.2f}": {
                "W_insertion_N": wedge_force(P_EB, g.alpha_deg, mu),
                "W_retention_N": wedge_force(P_EB, g.beta_deg, mu),
            }
            for mu in (0.0, 0.1, 0.2, 0.3)
        },
    }


def main() -> None:
    g, m = Geometry(), Material()

    print("=" * 74)
    print("ANALYTICAL REFERENCE")
    print("=" * 74)
    print(f"  Units            : mm - tonne - s - N - MPa")
    print(f"  Material         : {m.name}")
    print(f"  kappa convention : {KAPPA_CONVENTION}")

    print("\n  SELF-CHECKS")
    for line in _self_check(g, m):
        print(f"    [ok] {line}")

    tgt = build_targets(g, m)
    s1, s2 = tgt["stage1"], tgt["stage2"]

    print("\n" + "-" * 74)
    print("  GEOMETRY")
    print("-" * 74)
    print(f"    L = {g.L:.2f} mm   b = {g.b:.2f} mm   t = {g.t:.2f} mm   L/t = {g.slenderness:.1f}")
    print(f"    y = {g.y:.2f} mm   (= tooth protrusion = notch depth)")
    print(f"    alpha = {g.alpha_deg:.0f} deg, beta = {g.beta_deg:.0f} deg"
          f"   [from the insertion axis]")
    print(f"    I = {s1['I_mm4']:.6f} mm^4")

    print("\n" + "-" * 74)
    print("  STAGE 1 - idealized root clamp, linear elastic, no contact")
    print("-" * 74)
    print(f"    Euler-Bernoulli           P_EB = {s1['P_euler_bernoulli_N']:.4f} N")
    print(f"    shear correction                 {s1['shear_correction_pct']:.4f} %")
    print(f"    Timoshenko                P_T  = {s1['P_timoshenko_N']:.4f} N   <-- FE comparison")
    print(f"    acceptance band (+/-2%)          "
          f"{s1['acceptance_band_N'][0]:.3f} .. {s1['acceptance_band_N'][1]:.3f} N")
    print(f"    root stress (linear)             {s1['root_stress_MPa']:.2f} MPa")
    print(f"    root strain                      {s1['root_strain_pct']:.4f} %")
    print(f"    permissible strain               {s1['permissible_strain_pct']:.2f} % "
          f"({m.permissible_strain_factor:.2f} x eps_break)")
    print(f"    utilisation                      {s1['strain_utilisation_pct']:.0f} %")
    print(f"    f1 = {s1['f1_Hz']:,.0f} Hz   T1 = {s1['T1_ms']:.4f} ms"
          f"   c = {s1['wave_speed_mm_s']:,.0f} mm/s")

    print("\n" + "-" * 74)
    print("  STAGE 1 - deflected shape, second check on the same model")
    print("-" * 74)
    print(f"    {'x [mm]':>8} {'x/L':>6} {'w_EB [mm]':>12} {'w_T [mm]':>12} {'diff':>8}")
    for key, st in s1["deflected_shape_mm"].items():
        w_eb, w_t = st["w_euler_bernoulli"], st["w_timoshenko"]
        diff = (w_t / w_eb - 1.0) * 100.0 if w_eb else 0.0
        print(f"    {st['x_over_L'] * g.L:>8.1f} {st['x_over_L']:>6.2f} "
              f"{w_eb:>12.4f} {w_t:>12.4f} {diff:>7.2f} %")
    print(f"\n    The tip deflection is prescribed in the FE model, so the shape between")
    print(f"    the root and the tip is an independent check: it tests the response")
    print(f"    without the reaction force entering the comparison. Shear moves the")
    print(f"    shape by about 1 % at quarter span, which is at the edge of what the")
    print(f"    mesh sequence can resolve - the shape check discriminates beam theory")
    print(f"    from a wrong model, not Timoshenko from Euler-Bernoulli.")

    print("\n" + "-" * 74)
    print("  STAGE 2 - wedge relation, idealized analytical reference")
    print("-" * 74)
    print(f"    {'mu':>6} {'W_insertion [N]':>17} {'W_retention [N]':>17}")
    for mu in (0.0, 0.1, 0.2, 0.3):
        k = f"mu_{mu:.2f}"
        print(f"    {mu:>6.2f} {s2[k]['W_insertion_N']:>17.4f} {s2[k]['W_retention_N']:>17.4f}")
    print(f"\n    At mu = 0 and beta = 45 deg, tan(beta) = 1 and W_retention returns P_EB.")
    print(f"    The frictionless case separates the contact formulation from the")
    print(f"    friction model in the Stage 2 comparison.")

    print("\n" + "-" * 74)
    print("  MATERIAL MODEL - linear-elastic assumption")
    print("-" * 74)
    print(f"    initial tangent modulus (ISO 527)      {m.E:,.0f} MPa")
    print(f"    secant modulus at break ({m.sigma_break:.0f}/{m.eps_break:.3f})    "
          f"{m.E_secant_at_break:,.0f} MPa")
    print(f"    difference over the range to break     "
          f"{(1 - m.E_secant_at_break / m.E) * 100:.0f} %")
    print(f"    linear stress at design strain         {s1['root_stress_MPa']:.1f} MPa "
          f"= {s1['root_stress_MPa'] / m.sigma_break * 100:.0f} % of break stress")
    print(f"    design strain                          "
          f"{s1['root_strain_pct'] / (m.eps_break * 100) * 100:.0f} % of break strain")
    print(f"\n    Tangent and secant moduli differ over the range to break, so the")
    print(f"    response is non-linear. The linear-elastic model should therefore be")
    print(f"    regarded as a verification model rather than a prediction of the")
    print(f"    physical snap-fit behaviour. The Stage 1 comparison is unaffected: a")
    print(f"    linear FE model is compared against a linear closed form, which tests")
    print(f"    the numerical implementation and not the material model.")

    print("\n" + "-" * 74)
    print("  CREEP / RELAXATION")
    print("-" * 74)
    print(f"    Datasheet: {m.E:,.0f} MPa tensile modulus, {m.E_creep_1000h:,.0f} MPa")
    print(f"    tensile creep modulus after 1000 h at 23 C (ISO 899-1).")
    print(f"    The lance is held at approximately constant strain, so stress")
    print(f"    relaxation is the applicable quantity; the creep modulus is an")
    print(f"    indicator only. No retention force after 1000 h is reported.")

    print("\n  Poisson sensitivity (nu assumed, not from the datasheet):")
    for n in (0.30, 0.35, 0.40):
        mt = Material(nu=n)
        print(f"    nu = {n:.2f}  ->  P_T = {tip_force_timoshenko(g, mt):.4f} N")
    print("    Spread below 0.1 %; the assumed value does not affect the reference.")

    out = Path(__file__).resolve().parent.parent / "results" / "analytical_targets.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(tgt, indent=2))
    print(f"\n  Targets written to: {out.relative_to(out.parents[1])}")
    print("=" * 74)


if __name__ == "__main__":
    main()
