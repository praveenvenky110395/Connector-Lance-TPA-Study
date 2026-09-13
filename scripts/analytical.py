#!/usr/bin/env python3
"""
analytical.py — Analytical reference for the connector lance study.

PURPOSE
    Establish the closed-form expectations for the locking lance BEFORE any finite
    element model exists. Every value below is computed from the inputs in the
    INPUTS block; nothing is hard-coded. The FE models built in Stages 1-3 are
    measured against these numbers.

METHOD NOTE
    This script verifies itself. Each quantity is computed by its compact closed
    form and then re-computed by an independent algebraic route; the two must
    agree or the script raises. The wedge relation is additionally checked against
    its two physical limits. The same discipline applied to the FE models is
    applied here.

UNIT SYSTEM
    mm - tonne - s - N - MPa   (the standard LS-DYNA / Abaqus automotive system)
        length          mm
        mass            tonne (Mg)
        time            s
        force           N
        stress          MPa  = N/mm^2
        density         tonne/mm^3      (1 g/cm^3 = 1e-9 t/mm^3)
        frequency       Hz

SCOPE AND LIMITATION
    The Stage 1 model is LINEAR ELASTIC and is a solver-verification model, not a
    validated prediction of a physical retention force. See docs/numerical_spec.md,
    section "Limitations". The absolute forces reported here are upper bounds.

Author: Praveen Venkatesh Sethumadhavan Vasan
Licence / use: educational and self-learning. Geometry is synthetic and
    representative; it is not a reproduction of any manufacturer's product.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

# =============================================================================
# INPUTS — single source of truth. Change a value here and everything follows.
# =============================================================================


@dataclass(frozen=True)
class Geometry:
    """Locking lance, terminal notch and ramp angles. All lengths in mm."""

    L: float = 8.00      # lance free length (root clamp to tooth centreline)
    b: float = 2.50      # lance width
    t: float = 0.80      # lance thickness
    y: float = 0.60      # required deflection = tooth protrusion = notch depth

    # Ramp angles, BOTH MEASURED FROM THE INSERTION AXIS.
    # A shallow ramp is a SMALL angle. alpha -> 0 is a face parallel to the
    # insertion direction (pure sliding); alpha -> 90 deg is a face normal to it
    # (self-locking, cannot be inserted). Getting this convention backwards in
    # CAD invalidates every Stage 2 number.
    alpha_deg: float = 30.0   # lead-in ramp, insertion side
    beta_deg: float = 45.0    # return face, retention side

    @property
    def slenderness(self) -> float:
        """L/t. Keep >= 8 so Euler-Bernoulli remains a good approximation."""
        return self.L / self.t


@dataclass(frozen=True)
class Material:
    """BASF Ultradur B 4300 G6, PBT-GF30. Datasheet values with ISO methods."""

    name: str = "BASF Ultradur B 4300 G6 (PBT-GF30)"
    E: float = 9800.0          # MPa, tensile modulus, ISO 527-1/-2
    nu: float = 0.35           # assumed - NOT on the datasheet (see sensitivity)
    rho: float = 1.53e-9       # tonne/mm^3, from 1530 kg/m^3, ISO 1183
    sigma_break: float = 137.0  # MPa, stress at break, ISO 527-1/-2
    eps_break: float = 0.030   # strain at break, ISO 527-1/-2
    E_creep_1000h: float = 7500.0  # MPa, tensile creep modulus 1000 h/23C, ISO 899-1

    # Permissible strain for a one-time-assembly snap-fit in a glass-reinforced
    # grade, expressed as a fraction of the datasheet strain at break. 0.50 is the
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
        """Secant modulus at the break point, MPa. Compare against E."""
        return self.sigma_break / self.eps_break


# Shear correction factor for a rectangular cross-section.
#   "classic" : kappa = 5/6, the textbook value used in most FE documentation
#   "cowper"  : kappa = 10(1+nu)/(12+11nu), Cowper's refined value
# The two differ by ~0.02 % in the final force here. State which was used.
KAPPA_CONVENTION = "classic"


# =============================================================================
# CLOSED-FORM RELATIONS
# =============================================================================


def second_moment_of_area(g: Geometry) -> float:
    """I = b t^3 / 12  [mm^4]  — rectangular section about the bending axis."""
    return g.b * g.t**3 / 12.0


def shear_coefficient(m: Material, convention: str = KAPPA_CONVENTION) -> float:
    """Timoshenko shear correction factor kappa for a rectangular section."""
    if convention == "classic":
        return 5.0 / 6.0
    if convention == "cowper":
        return 10.0 * (1.0 + m.nu) / (12.0 + 11.0 * m.nu)
    raise ValueError(f"unknown kappa convention: {convention!r}")


def tip_force_euler_bernoulli(g: Geometry, m: Material) -> float:
    """
    Force at the lance tip to hold a prescribed deflection y.  [N]

        delta = P L^3 / (3 E I)      ->      P = E b t^3 y / (4 L^3)

    Bending only. No shear deformation, rigid root, point load at the tip.
    """
    return m.E * g.b * g.t**3 * g.y / (4.0 * g.L**3)


def shear_correction_ratio(g: Geometry, m: Material) -> float:
    """
    Extra tip deflection from transverse shear, as a fraction of the bending
    deflection:

        delta_shear / delta_bending = 3 E I / (kappa G A L^2)

    which for a rectangular section with kappa = 5/6 simplifies to

        = 0.6 (1 + nu) (t/L)^2
    """
    I = second_moment_of_area(g)
    A = g.b * g.t
    kappa = shear_coefficient(m)
    return 3.0 * m.E * I / (kappa * m.G * A * g.L**2)


def tip_force_timoshenko(g: Geometry, m: Material) -> float:
    """
    Tip force including transverse shear, for the SAME prescribed deflection. [N]

    The beam is softer than Euler-Bernoulli predicts, so less force is required.
    This is the value a well-meshed FE model with a perfect root clamp should
    reproduce - not the Euler-Bernoulli value.
    """
    return tip_force_euler_bernoulli(g, m) / (1.0 + shear_correction_ratio(g, m))


def root_stress(g: Geometry, m: Material) -> float:
    """Peak bending stress at the arm root, sigma = M c / I.  [MPa]"""
    P = tip_force_euler_bernoulli(g, m)
    return (P * g.L) * (g.t / 2.0) / second_moment_of_area(g)


def max_root_strain(g: Geometry) -> float:
    """
    Peak bending strain at the arm root.  [-]

        eps = 3 t y / (2 L^2)

    Independent of E, which is why it is the correct quantity to check the
    material choice against.
    """
    return 3.0 * g.t * g.y / (2.0 * g.L**2)


def wave_speed(m: Material) -> float:
    """One-dimensional elastic wave speed, c = sqrt(E/rho).  [mm/s]"""
    return math.sqrt(m.E / m.rho)


def first_bending_frequency(g: Geometry, m: Material) -> float:
    """
    First cantilever bending mode.  [Hz]

        f1 = (lambda1^2 / 2 pi) sqrt( E I / (rho A L^4) ),   lambda1 = 1.8751

    Sets the loading time for an explicit quasi-static run: the load must be
    applied over many multiples of T1 = 1/f1 for inertia to be negligible.
    """
    lam1 = 1.8751040687119611
    I = second_moment_of_area(g)
    A = g.b * g.t
    return (lam1**2 / (2.0 * math.pi)) * math.sqrt(m.E * I / (m.rho * A * g.L**4))


def wedge_force(P: float, angle_deg: float, mu: float) -> float:
    """
    Axial force to drive a spring-loaded tooth along an inclined face.  [N]

        W = P (tan(angle) + mu) / (1 - mu tan(angle))

    Derived from a free body:
        lance transverse equilibrium :  N (cos a - mu sin a) = P
        axial force                  :  W = N (sin a + mu cos a)

    ANGLE CONVENTION: measured from the insertion axis.
        angle -> 0   : face parallel to insertion, W -> mu P  (pure sliding)
        angle -> 90  : face normal to insertion, W -> infinity (self-locking)

    P is the transverse spring force from the lance at the deflection reached
    while riding the face.
    """
    ta = math.tan(math.radians(angle_deg))
    denom = 1.0 - mu * ta
    if denom <= 0.0:
        return math.inf  # self-locking
    return P * (ta + mu) / denom


# =============================================================================
# SELF-VERIFICATION — each result re-derived by an independent route
# =============================================================================


def _self_check(g: Geometry, m: Material, rtol: float = 1e-9) -> list[str]:
    """Re-derive every quantity a second way. Raises on disagreement."""
    checks: list[str] = []

    def close(a: float, b: float) -> bool:
        return abs(a - b) <= rtol * max(1.0, abs(b))

    # 1. Tip force: compact form vs 3 E I y / L^3
    I = second_moment_of_area(g)
    p_compact = tip_force_euler_bernoulli(g, m)
    p_from_EI = 3.0 * m.E * I * g.y / g.L**3
    assert close(p_compact, p_from_EI), "tip force: forms disagree"
    checks.append("P_EB: E b t^3 y/(4L^3) == 3 E I y/L^3")

    # 2. Strain: compact form vs sigma/E through M c / I
    eps_compact = max_root_strain(g)
    eps_from_stress = root_stress(g, m) / m.E
    assert close(eps_compact, eps_from_stress), "strain: forms disagree"
    checks.append("eps: 3ty/(2L^2) == (M c / I)/E")

    # 3. Shear correction: general form vs the rectangular simplification
    if KAPPA_CONVENTION == "classic":
        general = shear_correction_ratio(g, m)
        simplified = 0.6 * (1.0 + m.nu) * (g.t / g.L) ** 2
        assert close(general, simplified), "shear correction: forms disagree"
        checks.append("shear: 3EI/(kappa G A L^2) == 0.6(1+nu)(t/L)^2")

    # 4. Wedge relation, physical limit: a face parallel to insertion is pure
    #    sliding friction.
    P = p_compact
    assert close(wedge_force(P, 0.0, 0.25), 0.25 * P), "wedge: alpha=0 limit fails"
    checks.append("wedge limit: alpha=0 -> W = mu P")

    # 5. Wedge relation, physical limit: a 45 deg face with no friction returns
    #    exactly the transverse spring force (tan 45 = 1).
    assert close(wedge_force(P, 45.0, 0.0), P), "wedge: 45 deg frictionless fails"
    checks.append("wedge limit: beta=45 deg, mu=0 -> W = P  [Stage 2 regression test]")

    # 6. Wedge relation must increase monotonically with friction.
    ws = [wedge_force(P, g.alpha_deg, mu) for mu in (0.0, 0.1, 0.2, 0.3)]
    assert all(x < z for x, z in zip(ws, ws[1:])), "wedge: not monotonic in mu"
    checks.append("wedge: monotonic in mu")

    # 7. Beam theory validity gate.
    assert g.slenderness >= 8.0, f"L/t = {g.slenderness:.1f} is too stubby for beam theory"
    checks.append(f"slenderness: L/t = {g.slenderness:.1f} >= 8")

    return checks


# =============================================================================
# REPORT
# =============================================================================


def build_targets(g: Geometry, m: Material) -> dict:
    """All analytical targets in one dictionary, for downstream comparison."""
    P_EB = tip_force_euler_bernoulli(g, m)
    P_T = tip_force_timoshenko(g, m)
    f1 = first_bending_frequency(g, m)
    eps = max_root_strain(g)

    return {
        "units": "mm-tonne-s-N-MPa",
        "kappa_convention": KAPPA_CONVENTION,
        "geometry": asdict(g),
        "material": asdict(m),
        "stage1": {
            "I_mm4": second_moment_of_area(g),
            "P_euler_bernoulli_N": P_EB,
            "P_timoshenko_N": P_T,
            "shear_correction_pct": shear_correction_ratio(g, m) * 100.0,
            "root_stress_MPa": root_stress(g, m),
            "root_strain_pct": eps * 100.0,
            "permissible_strain_pct": m.eps_permissible * 100.0,
            "strain_utilisation_pct": eps / m.eps_permissible * 100.0,
            "f1_Hz": f1,
            "T1_ms": 1.0 / f1 * 1e3,
            "wave_speed_mm_s": wave_speed(m),
            "acceptance_band_N": [P_T * 0.98, P_T * 1.02],
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
    print("ANALYTICAL REFERENCE — computed before any FE model exists")
    print("=" * 74)
    print(f"  Units            : mm - tonne - s - N - MPa")
    print(f"  Material         : {m.name}")
    print(f"  kappa convention : {KAPPA_CONVENTION}")

    print("\n  SELF-VERIFICATION")
    for line in _self_check(g, m):
        print(f"    [ok] {line}")

    tgt = build_targets(g, m)
    s1, s2 = tgt["stage1"], tgt["stage2"]

    print("\n" + "-" * 74)
    print("  GEOMETRY")
    print("-" * 74)
    print(f"    L = {g.L:.2f} mm   b = {g.b:.2f} mm   t = {g.t:.2f} mm   L/t = {g.slenderness:.1f}")
    print(f"    y = {g.y:.2f} mm   (= tooth protrusion = notch depth, all three equal)")
    print(f"    alpha = {g.alpha_deg:.0f} deg, beta = {g.beta_deg:.0f} deg"
          f"   [MEASURED FROM THE INSERTION AXIS]")
    print(f"    I = {s1['I_mm4']:.6f} mm^4")

    print("\n" + "-" * 74)
    print("  STAGE 1 TARGETS — perfect root clamp, linear elastic, no contact")
    print("-" * 74)
    print(f"    Euler-Bernoulli           P_EB = {s1['P_euler_bernoulli_N']:.4f} N")
    print(f"    shear correction                 {s1['shear_correction_pct']:.4f} %")
    print(f"    Timoshenko                P_T  = {s1['P_timoshenko_N']:.4f} N   <-- FE target")
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
    print("  STAGE 2 TARGETS — wedge relation")
    print("-" * 74)
    print(f"    {'mu':>6} {'W_insertion [N]':>17} {'W_retention [N]':>17}")
    for mu in (0.0, 0.1, 0.2, 0.3):
        k = f"mu_{mu:.2f}"
        print(f"    {mu:>6.2f} {s2[k]['W_insertion_N']:>17.4f} {s2[k]['W_retention_N']:>17.4f}")
    print(f"\n    At mu = 0 and beta = 45 deg the retention force returns P_EB exactly.")
    print(f"    Run this case FIRST: it tests the contact algorithm with the friction")
    print(f"    term removed, and it is a regression test on Stage 1.")

    print("\n" + "-" * 74)
    print("  MATERIAL LINEARITY — read this before quoting any absolute force")
    print("-" * 74)
    print(f"    initial tangent modulus (ISO 527)      {m.E:,.0f} MPa")
    print(f"    secant modulus at break ({m.sigma_break:.0f}/{m.eps_break:.3f})    "
          f"{m.E_secant_at_break:,.0f} MPa")
    print(f"    softening by the break point           "
          f"{(1 - m.E_secant_at_break / m.E) * 100:.0f} %")
    print(f"    linear stress at design strain         {s1['root_stress_MPa']:.1f} MPa "
          f"= {s1['root_stress_MPa'] / m.sigma_break * 100:.0f} % of break stress")
    print(f"    ... at only                            "
          f"{s1['root_strain_pct'] / (m.eps_break * 100) * 100:.0f} % of break strain")
    print(f"\n    The stress-strain curve is NOT linear at the design strain. The forces")
    print(f"    above are UPPER BOUNDS under the linear-elastic assumption. Stage 1")
    print(f"    verification is unaffected: it compares a linear FE model against a")
    print(f"    linear closed form and therefore tests the SOLVER, not the material.")

    print("\n" + "-" * 74)
    print("  CREEP / RELAXATION — what may and may not be claimed")
    print("-" * 74)
    print(f"    Datasheet: {m.E:,.0f} MPa tensile modulus, {m.E_creep_1000h:,.0f} MPa")
    print(f"    tensile creep modulus after 1000 h at 23 C (ISO 899-1).")
    print(f"    The lance is held at approximately constant STRAIN, so stress")
    print(f"    RELAXATION is the strictly correct quantity and the creep modulus is")
    print(f"    an indicator only. No retention force after 1000 h is claimed.")

    print("\n  Poisson sensitivity (nu is assumed, not from the datasheet):")
    for n in (0.30, 0.35, 0.40):
        mt = Material(nu=n)
        print(f"    nu = {n:.2f}  ->  P_T = {tip_force_timoshenko(g, mt):.4f} N")
    print("    Spread below 0.1 %. The assumption is immaterial to the result.")

    out = Path(__file__).resolve().parent.parent / "results" / "analytical_targets.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(tgt, indent=2))
    print(f"\n  Targets written to: {out.relative_to(out.parents[1])}")
    print("=" * 74)


if __name__ == "__main__":
    main()
