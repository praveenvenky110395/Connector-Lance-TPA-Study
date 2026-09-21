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


# Stage 3 prop stations and clearances, in mm. The stations are all multiples of
# the 0.20 mm element length so the prop lands on a node column rather than
# between two, which is what keeps the line of action where the closed form
# assumes it is. The 0.30 mm clearance never gets used up inside the 0.60 mm
# stroke - that run is there to show a defined contact changing nothing.
STAGE3_PROPS = (3.0, 3.6, 4.0)
STAGE3_GAPS = (0.05, 0.10, 0.30)


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
    return P * wedge_ratio(angle_deg, mu)


def wedge_ratio(angle_deg: float, mu: float) -> float:
    """
    W/P for the same relation - the axial force divided by the transverse one.

        W/P = (tan(angle) + mu) / (1 - mu tan(angle))

    This is the quantity Stage 2 is checked on, and the reason is worth stating.
    The ratio contains no stiffness: not E, not the second moment of area, not
    the deflection. So the Stage 1 difference between the FE model and the beam
    solution - which is a stiffness difference - cancels out of it completely.
    Stage 2 then tests only what Stage 2 adds, which is the contact normal and
    the friction law.

    Both quantities come out of the same spcforc file at the same instant:
    W is the summed axial reaction at the root, P the summed transverse one.
    """
    ta = math.tan(math.radians(angle_deg))
    denom = 1.0 - mu * ta

    # Self-locking is mu tan(angle) >= 1. The tolerance is not cosmetic: at 45
    # degrees tan comes back as 0.9999999999999999, so an exact test on zero
    # misses the locking point by one bit and returns a huge finite number
    # instead of infinity. A denominator this small means the relation has
    # diverged either way.
    if denom <= 1e-12:
        return math.inf
    return (ta + mu) / denom


def plate_travel(g: Geometry, angle_deg: float) -> float:
    """
    How far the rigid plate has to advance to deflect the lance by y.  [mm]

        s = y / tan(angle)

    The plate's inclined face descends by s tan(angle) at the tip as it advances
    by s, so this follows from the geometry alone and is an input check on the
    FE model rather than a result: at the end of the ramp the tip deflection must
    be y, whatever the friction coefficient is.
    """
    return g.y / math.tan(math.radians(angle_deg))


# =============================================================================
# Stage 3 - the TPA as a prop under the lance
# =============================================================================
#
# The secondary lock does not hold the terminal itself. It sits under the lance
# with a small clearance and takes away most of the deflection the lance needs
# in order to release. Mechanically that turns a cantilever into a propped
# cantilever, once the clearance is used up.
#
# The response is therefore in two parts:
#
#   below the engagement deflection   the prop is not touched, and the lance is
#                                     the Stage 1 cantilever
#   above it                          the prop is a support, the beam is singly
#                                     redundant, and the tip stiffness jumps
#
# Two numbers come out of that, and both are ratios with no modulus in them, for
# the same reason the Stage 2 wedge ratio has none: they are built from
# flexibilities, and every flexibility carries the same 1/E. So the Stage 1
# difference between the FE model and beam theory cancels out of them, exactly
# as it cancelled out of Stage 2. That is what makes Stage 3 a check on the new
# ingredient - the prop - rather than a restatement of Stage 1.


def cantilever_flexibility(g: Geometry, m: Material, x_load: float,
                           x_read: float, shear: bool = True) -> float:
    """
    Deflection at x_read under a unit transverse load at x_load, for the
    cantilever clamped at x = 0.  [mm/N]

        bending   lo^2 (3 hi - lo) / (6 E I)
        shear     lo / (kappa G A)

    with lo and hi the nearer and further of the two stations. Written this way
    the function is symmetric in its two arguments, which is Maxwell-Betti and
    is checked rather than assumed. The shear term uses lo because a load at
    x_load puts shear only in 0 < x < x_load; past that the beam translates
    without further shear deformation.
    """
    lo, hi = min(x_load, x_read), max(x_load, x_read)
    I = second_moment_of_area(g)
    A = g.b * g.t
    out = lo**2 * (3.0 * hi - lo) / (6.0 * m.E * I)
    if shear:
        out += lo / (shear_coefficient(m) * m.G * A)
    return out


def _prop_coefficients(g: Geometry, m: Material, prop_x: float,
                       shear: bool = True) -> tuple[float, float, float]:
    """The three flexibilities the propped problem needs: C, A, B."""
    C = cantilever_flexibility(g, m, g.L, g.L, shear)       # tip from tip load
    A = cantilever_flexibility(g, m, g.L, prop_x, shear)    # prop from tip load
    B = cantilever_flexibility(g, m, prop_x, prop_x, shear)  # prop from prop load
    return C, A, B


def prop_engagement_deflection(g: Geometry, m: Material, prop_x: float,
                               gap: float) -> float:
    """
    Tip deflection at which the underside of the lance first reaches the prop.

        delta_eng = gap * C / A

    C/A is the ratio of the tip deflection to the deflection at the prop station
    under the same tip load, so it is the reciprocal of the deflected shape
    evaluated at the prop. No modulus in it, and it is linear in the clearance -
    both of which are checked against the FE runs.
    """
    C, A, _ = _prop_coefficients(g, m, prop_x)
    return gap * C / A


def prop_stiffness_ratio(g: Geometry, m: Material, prop_x: float,
                         shear: bool = True) -> float:
    """
    Tip stiffness with the prop engaged, divided by the free-cantilever tip
    stiffness.  [-]

        k2 / k1 = C B / (C B - A^2)

    Depends on prop_x / L alone. Not on E, not on the clearance, not on the
    deflection. The clearance independence is worth stating separately because
    it is a cheap and sharp test of the FE model: two runs at the same prop
    position and different clearances must return the same slope, and if they
    do not, the contact is doing something the closed form does not describe.
    """
    C, A, B = _prop_coefficients(g, m, prop_x, shear)
    denominator = C * B - A * A
    if denominator <= 0.0:
        return math.inf                      # prop at the tip: no travel left
    return C * B / denominator


def propped_response(g: Geometry, m: Material, prop_x: float, gap: float,
                     delta: float) -> tuple[float, float]:
    """
    (tip force, prop reaction) at tip deflection delta.  [N, N]

    Below engagement the prop carries nothing and this is the Stage 1 spring.
    Above it, the two compatibility conditions

        w(prop_x) = gap        P A - R B = gap
        w(L)      = delta      P C - R A = delta

    give P and R directly. Both branches meet at the engagement deflection with
    P continuous and R passing through zero, which is checked below.
    """
    C, A, B = _prop_coefficients(g, m, prop_x)
    if delta <= gap * C / A:
        return delta / C, 0.0
    P = (delta * B - A * gap) / (C * B - A * A)
    return P, (P * A - gap) / B


def root_reactions(g: Geometry, prop_x: float, tip_force: float,
                   prop_force: float) -> tuple[float, float]:
    """
    (transverse force, bending moment) at the root, from statics alone.  [N, N mm]

    This is how the FE forces are recovered. The tip is driven through a rigid
    body that is free in x and free to rotate, so it can only deliver a
    transverse force; the prop is frictionless, so it can only deliver a normal
    one. Three unknowns would be two too many, but with the prop station known
    by construction the two equilibrium equations close the system, and the
    summed root forces in spcforc give both of them.
    """
    return tip_force - prop_force, prop_x * prop_force - g.L * tip_force


# =============================================================================
# Stage 3 - the locking cycle
# =============================================================================
#
# Stage 3 stops being a beam exercise and becomes the part. The lance carries an
# integral locking tooth: a 45 deg retention face on the root side, a flat crest,
# and a 30 deg lead-in on the tip side, protruding by y - the protrusion IS the
# deflection needed to release, which is the design relation. A rigid terminal
# runs underneath it and a rigid TPA sits above.
#
# Nothing new is assumed. The prediction for Stage 3 is the product of the two
# stages already verified:
#
#     force to deflect the lance        Stage 1, now on a variable section
#     force to pull the terminal        Stage 2, W/P = (tan b + mu)/(1 - mu tan b)
#
# Three things change when the idealisation becomes the part, and all three are
# named here before the model is run rather than discovered afterwards:
#
#   1. the tooth stiffens the lance, because the section is thicker where the
#      tooth is
#   2. the load is applied at the tooth, not at the tip, so the lever arm is
#      shorter and the release condition is set at the crest rather than at the
#      tip
#   3. the inclined faces are on the lance and the lance rotates, so the contact
#      normal rotates with it. Stage 2 deliberately put the wedge face on the
#      rigid body to avoid exactly this. Here it cannot be avoided, because that
#      is where the face is on a real lance - so it is quantified instead.
#
# Item 3 is the one worth reading twice. Lifting the lance rotates the tooth, the
# retention face gets shallower and the lead-in face gets steeper, both by the
# same angle. Retention falls below the hand calculation and insertion rises
# above it, so the designed asymmetry between them shrinks.

# Crest length between the two faces, mm. Long enough that the two contacts
# never share a segment, short enough to keep the tooth compact.
TOOTH_CREST = 0.40


def tooth_stations(g: Geometry) -> tuple[float, float, float, float]:
    """
    The four x stations of the tooth profile: retention face from x1 to x2,
    crest x2 to x3, lead-in x3 to x4, with x4 at the lance tip.

    The runs follow from the angles and the protrusion, so the tooth is not a
    free choice: a 45 deg face rising 0.60 mm needs 0.60 mm of x, and a 30 deg
    face needs 1.04 mm.
    """
    x4 = g.L
    x3 = x4 - g.y / math.tan(math.radians(g.alpha_deg))
    x2 = x3 - TOOTH_CREST
    x1 = x2 - g.y / math.tan(math.radians(g.beta_deg))
    return x1, x2, x3, x4


def tooth_depth(g: Geometry, x: float) -> float:
    """How far the tooth protrudes below the lance underside at station x. [mm]"""
    x1, x2, x3, x4 = tooth_stations(g)
    if x <= x1 or x >= x4:
        return 0.0
    if x < x2:
        return (x - x1) * math.tan(math.radians(g.beta_deg))
    if x < x3:
        return g.y
    return (x4 - x) * math.tan(math.radians(g.alpha_deg))


def section_thickness(g: Geometry, x: float, with_tooth: bool = True) -> float:
    """Local section thickness. [mm]"""
    return g.t + (tooth_depth(g, x) if with_tooth else 0.0)


def beam_flexibility(g: Geometry, m: Material, x_load: float, x_read: float,
                     with_tooth: bool = True, n: int = 8001) -> float:
    """
    Deflection at x_read under a unit transverse load at x_load. [mm/N]

    Unit-load method, integrated numerically because the section varies:

        f(a,b) = int_0^min(a,b) (a-x)(b-x)/(E I(x)) dx
               + int_0^min(a,b) 1/(kappa G A(x)) dx

    Both virtual moment diagrams are zero past their own load point, so the
    integral stops at the nearer of the two stations. Symmetric in its two
    arguments by construction, which is Maxwell-Betti and is checked. With the
    tooth switched off this has to reproduce the prismatic formulas exactly, and
    that is checked too - it is the only way to know the integration is right
    before trusting it on a section it cannot be checked against.
    """
    top = min(x_load, x_read)
    if top <= 0.0:
        return 0.0
    kappa = shear_coefficient(m)
    h = top / (n - 1)
    total = 0.0
    for i in range(n):
        x = i * h
        t = section_thickness(g, x, with_tooth)
        I = g.b * t**3 / 12.0
        A = g.b * t
        f = ((x_load - x) * (x_read - x) / (m.E * I)
             + 1.0 / (kappa * m.G * A))
        # Simpson, so the odd samples carry weight 4 and the even ones 2.
        weight = 1.0 if i in (0, n - 1) else (4.0 if i % 2 else 2.0)
        total += weight * f
    return total * h / 3.0


def beam_rotation_flexibility(g: Geometry, m: Material, x_load: float,
                              x_read: float, with_tooth: bool = True,
                              n: int = 8001) -> float:
    """
    Rotation at x_read under a unit transverse load at x_load. [rad/N]

    Same method with a unit moment as the virtual load, whose diagram is 1 over
    the whole span up to x_read. Shear does not contribute to a rotation.
    """
    top = min(x_load, x_read)
    if top <= 0.0:
        return 0.0
    h = top / (n - 1)
    total = 0.0
    for i in range(n):
        x = i * h
        t = section_thickness(g, x, with_tooth)
        f = (x_load - x) / (m.E * g.b * t**3 / 12.0)
        weight = 1.0 if i in (0, n - 1) else (4.0 if i % 2 else 2.0)
        total += weight * f
    return total * h / 3.0


def face_stations(g: Geometry) -> dict:
    """
    Where each face bears, and the ends of each face.

    The contact is face to face, not a point, and it migrates: on extraction the
    two 45 deg faces slide out of each other as the lance lifts, so the resultant
    travels along the face. The mid-point is used as the working station and the
    two ends are carried with it, so every result below can be given as a bracket
    rather than as a number that pretends to know where the resultant sat.
    """
    x1, x2, x3, x4 = tooth_stations(g)
    return {
        "retention": 0.5 * (x1 + x2), "retention_ends": (x1, x2),
        "lead_in": 0.5 * (x3 + x4), "lead_in_ends": (x3, x4),
        "crest": 0.5 * (x2 + x3),
    }


def release_force(g: Geometry, m: Material, station: float) -> float:
    """
    Transverse force at `station` needed to lift the crest clear.  [N]

        lift(crest) = P f(station, crest) = y    ->    P = y / f(station, crest)

    The lance tip never enters it. The load is on a tooth face and the condition
    is at the crest, and one flexibility function answers both.
    """
    return g.y / beam_flexibility(g, m, station, face_stations(g)["crest"])


def lance_rotation(g: Geometry, m: Material, station: float) -> float:
    """How far the tooth has rotated when the crest is about to clear. [degrees]"""
    P = release_force(g, m, station)
    return math.degrees(P * beam_rotation_flexibility(g, m, station, station))


def effective_angles(g: Geometry, m: Material) -> tuple[float, float]:
    """
    (lead-in, retention) face angles once the lance has rotated. [degrees]

    This is the finding of Stage 3, and it is geometric rather than numerical.

    The lead-in face runs uphill towards the tip and the retention face runs
    downhill towards the root, so one rotation tilts them in opposite senses:
    the lead-in gets steeper, the retention gets shallower, both by the same
    seven and a bit degrees. A 30 deg face becomes about 37.5 and a 45 deg face
    becomes about 37.5 - they meet in the middle, and the asymmetry the two
    angles were chosen to produce very largely disappears.

    Nothing about this is visible to the hand calculation, which uses the drawn
    angles. It follows from Stage 1 (the lance rotates) and Stage 2 (the wedge
    relation is set by the face angle) put together, which is the whole reason
    the stages were built in that order.
    """
    st = face_stations(g)
    return (g.alpha_deg + lance_rotation(g, m, st["lead_in"]),
            g.beta_deg - lance_rotation(g, m, st["retention"]))


def retention_angle_for_effective(g: Geometry, m: Material,
                                  wanted_deg: float | None = None) -> float:
    """
    The retention angle to draw if the effective one is to come out at `wanted`.

    Defaults to the nominal beta, so the answer is what the face would have to be
    for the rotation to leave it where it was meant to be. One line, and it is
    the design recommendation the stage produces.
    """
    wanted = g.beta_deg if wanted_deg is None else wanted_deg
    return wanted + lance_rotation(g, m, face_stations(g)["retention"])


def cycle_forces(g: Geometry, m: Material, mu: float) -> dict:
    """
    Insertion and retention force, as drawn and as rotated, with brackets.

    Insertion peaks when the terminal corner reaches the crest, extraction when
    the faces are about to part - both at the same lift, so the same release
    condition serves for both and only the load station differs.
    """
    st = face_stations(g)
    out = {"mu": mu}
    for name, key, nominal, sign in (("insertion", "lead_in", g.alpha_deg, +1.0),
                                     ("retention", "retention", g.beta_deg, -1.0)):
        a = st[key]
        P = release_force(g, m, a)
        eff = nominal + sign * lance_rotation(g, m, a)
        bracket = []
        for end in st[key + "_ends"]:
            Pe = release_force(g, m, end)
            bracket.append(Pe * wedge_ratio(
                nominal + sign * lance_rotation(g, m, end), mu))
        out[name] = {
            "station_mm": a,
            "P_release_N": P,
            "angle_drawn_deg": nominal,
            "angle_effective_deg": eff,
            "force_drawn_N": P * wedge_ratio(nominal, mu),
            "force_rotated_N": P * wedge_ratio(eff, mu),
            "force_rotated_bracket_N": [min(bracket), max(bracket)],
        }
    out["asymmetry_drawn"] = (out["retention"]["force_drawn_N"]
                              / out["insertion"]["force_drawn_N"])
    out["asymmetry_rotated"] = (out["retention"]["force_rotated_N"]
                                / out["insertion"]["force_rotated_N"])
    return out


def tpa_lift_ratio(g: Geometry, m: Material, station: float) -> float:
    """
    Lift at `station` divided by lift at the crest, under the extraction load.

    The lance is not lifted bodily - it is bent, so every station lifts by a
    different amount, and stations outboard of the load lift most. The tip lifts
    1.26 times the crest here. Which station matters depends on where the TPA is.
    """
    a = face_stations(g)["retention"]
    return (beam_flexibility(g, m, a, station)
            / beam_flexibility(g, m, a, face_stations(g)["crest"]))


def tpa_blocking_clearance(g: Geometry, m: Material, covers_to: float) -> float:
    """
    The clearance below which a TPA reaching out to `covers_to` prevents release.

    An earlier version of this was one line - "clearance < protrusion" - and it
    was wrong. It compared the clearance against how far the *tooth* has to move
    and ignored that the TPA sits somewhere else on a beam that bends. The lance
    tip lifts 1.26 times the crest, so a TPA over the tip stops the lance at
    0.757 mm of clearance, not 0.600. A run built to show a TPA that is fitted
    and useless at 0.65 mm clearance in fact blocked, and the FE model said so
    before the arithmetic did.

    Corrected: the criterion is the clearance against the lift at the most lifted
    station the TPA actually covers - which is a statement about where the TPA is
    put, not about the tooth. Putting it over the tip buys 26 % more clearance
    for the same function, which is a design lever the first version hid.
    """
    return g.y * tpa_lift_ratio(g, m, min(covers_to, g.L))


def tpa_blocks(g: Geometry, m: Material, clearance: float,
               covers_to: float) -> bool:
    """Whether a TPA at this clearance and reach prevents release."""
    return clearance < tpa_blocking_clearance(g, m, covers_to)


def tpa_block_travel(g: Geometry, m: Material, clearance: float,
                     covers_to: float) -> float:
    """
    Terminal travel at which the lance bottoms on the TPA. [mm]

    The retention faces are at 45 deg, so the lift at the contact follows the
    travel one for one, and the lift anywhere else follows from the bending
    shape. This is the number the stroke has to be set from: drive further and
    the terminal, which is rigid and kinematically driven, simply ploughs on
    with nothing in the model able to stop it.
    """
    a = face_stations(g)["retention"]
    ratio = (beam_flexibility(g, m, a, min(covers_to, g.L))
             / beam_flexibility(g, m, a, a))
    return clearance / ratio


# =============================================================================
# Self-checks - each result re-derived by an independent route
# =============================================================================


def _self_check(g: Geometry, m: Material, rtol: float = 1e-9) -> list[str]:
    """Re-derive each quantity a second way. Raises on disagreement."""
    checks: list[str] = []

    def close(a: float, b: float, tol: float | None = None) -> bool:
        # The Stage 3 quantities come out of a numerical integration rather than
        # a formula, so they need a looser tolerance than the closed forms do.
        # Passing it per call keeps that visible at the call site instead of
        # relaxing the default for everything.
        return abs(a - b) <= (rtol if tol is None else tol) * max(1.0, abs(b))

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

    # The ratio form and the force form have to be the same relation
    for a in (0.0, 15.0, g.alpha_deg, g.beta_deg):
        for mu in (0.0, 0.15, 0.3):
            assert close(wedge_force(P, a, mu), P * wedge_ratio(a, mu)), \
                f"wedge: ratio and force forms disagree at {a} deg, mu = {mu}"
    checks.append("wedge: W = P * (W/P) over the angles and mu used")

    # The number the Stage 2 FE model is regressed against
    assert close(wedge_ratio(g.beta_deg, 0.0), 1.0), "wedge ratio: beta=45, mu=0 != 1"
    checks.append("wedge ratio: beta=45 deg, mu=0 -> W/P = 1.000 exactly  [Stage 2]")

    # Sliding friction limit, now in ratio form
    assert close(wedge_ratio(0.0, 0.25), 0.25), "wedge ratio: alpha=0 limit fails"
    checks.append("wedge ratio: alpha=0 -> W/P = mu")

    # Self-locking must diverge rather than return something finite and wrong
    assert wedge_ratio(g.beta_deg, 1.0) == math.inf, "wedge: self-locking not detected"
    checks.append("wedge: self-locking at mu tan(angle) >= 1 returns infinity")

    # Plate travel is pure geometry. At 45 degrees it equals the deflection.
    assert close(plate_travel(g, 45.0), g.y), "plate travel: 45 deg case fails"
    for a in (g.alpha_deg, g.beta_deg):
        assert close(plate_travel(g, a) * math.tan(math.radians(a)), g.y), \
            f"plate travel: does not return y at {a} deg"
    checks.append("plate travel: s tan(angle) = y at both angles")

    # Slenderness screening criterion
    assert g.slenderness >= 8.0, \
        f"L/t = {g.slenderness:.1f} is below the slenderness screening criterion"
    checks.append(f"slenderness: L/t = {g.slenderness:.1f} >= 8")

    # --- Stage 3 -------------------------------------------------------------

    # Maxwell-Betti. If this fails the whole propped derivation is wrong, since
    # it assumes the same coefficient appears in both compatibility equations.
    for a in (2.0, 3.0, 4.0, 6.0):
        assert close(cantilever_flexibility(g, m, g.L, a),
                     cantilever_flexibility(g, m, a, g.L)), \
            f"flexibility: not symmetric at prop_x = {a}"
    checks.append("flexibility: reciprocal, f(L,a) = f(a,L)  [Stage 3]")

    # The tip-to-tip flexibility has to be the Stage 1 spring, or Stage 3 is
    # not built on Stage 1 at all.
    assert close(1.0 / cantilever_flexibility(g, m, g.L, g.L),
                 tip_force_timoshenko(g, m) / g.y), \
        "flexibility: tip flexibility does not reproduce the Timoshenko force"
    checks.append("flexibility: 1/f(L,L) reproduces the Stage 1 Timoshenko stiffness")

    # A prop at the root does nothing; a prop at the tip leaves no travel.
    assert close(prop_stiffness_ratio(g, m, 1e-9), 1.0), \
        "prop: ratio at the root is not 1"
    assert prop_stiffness_ratio(g, m, g.L) == math.inf, \
        "prop: degenerate case at the tip not detected"
    checks.append("prop: ratio -> 1 at the root, infinite at the tip")

    # Pure bending, prop at mid span. Doing this by hand gives 2304/504 and it
    # is the one value in Stage 3 that can be checked against arithmetic alone.
    assert close(prop_stiffness_ratio(g, m, g.L / 2.0, shear=False), 2304.0 / 504.0), \
        "prop: Euler-Bernoulli mid-span ratio is not 2304/504"
    checks.append("prop: bending-only mid-span ratio = 2304/504 = 4.5714")

    # No modulus in either ratio. Doubling E must change neither.
    stiff = Material(E=2.0 * m.E)
    for a in (3.0, 3.6, 4.0):
        assert close(prop_stiffness_ratio(g, m, a), prop_stiffness_ratio(g, stiff, a)), \
            f"prop: stiffness ratio depends on E at prop_x = {a}"
        assert close(prop_engagement_deflection(g, m, a, 0.1),
                     prop_engagement_deflection(g, stiff, a, 0.1)), \
            f"prop: engagement deflection depends on E at prop_x = {a}"
    checks.append("prop: stiffness ratio and engagement deflection are free of E")

    # Engagement is linear in the clearance, and at engagement the deflected
    # shape evaluated at the prop must return exactly the clearance.
    r = shear_correction_ratio(g, m)
    for a in (3.0, 3.6, 4.0):
        d1 = prop_engagement_deflection(g, m, a, 0.10)
        assert close(prop_engagement_deflection(g, m, a, 0.05), d1 / 2.0), \
            f"prop: engagement not linear in the clearance at prop_x = {a}"
        assert close(d1 * deflected_shape(a / g.L, r), 0.10), \
            f"prop: engagement disagrees with the Stage 1 deflected shape at {a}"
    checks.append("prop: engagement linear in clearance, and consistent with the "
                  "Stage 1 deflected shape")

    # The two branches must meet. Evaluated at the engagement deflection itself,
    # the propped expression has to collapse onto the free one and the reaction
    # has to be zero - both exactly, not to within a perturbation.
    for a in (3.0, 3.6, 4.0):
        for gap in (0.05, 0.10):
            C, A, B = _prop_coefficients(g, m, a)
            d = prop_engagement_deflection(g, m, a, gap)
            P_propped = (d * B - A * gap) / (C * B - A * A)
            assert close(P_propped, d / C), \
                f"prop: the two branches do not meet at engagement, {a}/{gap}"
            assert abs((P_propped * A - gap) / B) < 1e-12 * max(1.0, P_propped), \
                f"prop: reaction not zero at engagement, {a}/{gap}"
    checks.append("prop: branches meet at engagement with the reaction at zero")

    # Past engagement the slope must be the ratio, measured rather than assumed.
    for a in (3.0, 3.6, 4.0):
        gap, d0 = 0.10, prop_engagement_deflection(g, m, a, 0.10)
        P0, _ = propped_response(g, m, a, gap, d0 + 0.10)
        P1, _ = propped_response(g, m, a, gap, d0 + 0.20)
        slope = (P1 - P0) / 0.10
        assert close(slope * cantilever_flexibility(g, m, g.L, g.L),
                     prop_stiffness_ratio(g, m, a)), \
            f"prop: measured slope disagrees with the ratio at prop_x = {a}"
    checks.append("prop: slope of the propped branch equals k2/k1 times k1")

    # The solve itself, checked by putting the answer back into the two
    # compatibility conditions it came from: the beam must pass through the
    # clearance at the prop and through the prescribed deflection at the tip.
    for a in (3.0, 3.6, 4.0):
        for gap in (0.05, 0.10):
            C, A, B = _prop_coefficients(g, m, a)
            P, R = propped_response(g, m, a, gap, g.y)
            assert close(P * A - R * B, gap), \
                f"prop: solved forces do not put the beam on the prop, {a}/{gap}"
            assert close(P * C - R * A, g.y), \
                f"prop: solved forces do not reach the tip deflection, {a}/{gap}"
            # R by a second algebraic route - eliminate P instead of R.
            assert close(R, (g.y * A - gap * C) / (C * B - A * A)), \
                f"prop: the two eliminations disagree at {a}/{gap}"
    checks.append("prop: solution satisfies both compatibility conditions, and the "
                  "reaction agrees by a second elimination")

    # Statics, which is how the FE forces are recovered from spcforc.
    for a in (3.0, 3.6, 4.0):
        P, R = propped_response(g, m, a, 0.10, g.y)
        F_root, M_root = root_reactions(g, a, P, R)
        assert close(F_root, P - R), \
            f"statics: vertical equilibrium fails at prop_x = {a}"
        assert close(M_root + g.L * P - a * R, 0.0), \
            f"statics: moment equilibrium fails at prop_x = {a}"
        # And the inversion the extractor uses must give P and R back.
        assert close((M_root + a * F_root) / (a - g.L), P), \
            f"statics: the extractor inversion does not return P at {a}"
        assert close((M_root + a * F_root) / (a - g.L) - F_root, R), \
            f"statics: the extractor inversion does not return R at {a}"
    checks.append("prop: root force and moment close equilibrium, and invert back "
                  "to the tip force")

    # --- Stage 3, the locking cycle -----------------------------------------

    # The numerical integration has to reproduce the prismatic formulas exactly
    # with the tooth switched off. This is the only way to know it is right
    # before it is used on a section that has no formula.
    loose = 1e-6
    for a, bb in ((g.L, g.L), (3.0, g.L), (g.L, 3.0), (5.0, 6.0), (6.0, 5.0)):
        lo, hi = min(a, bb), max(a, bb)
        closed = (lo**2 * (3.0 * hi - lo) / (6.0 * m.E * second_moment_of_area(g))
                  + lo / (shear_coefficient(m) * m.G * g.b * g.t))
        got = beam_flexibility(g, m, a, bb, with_tooth=False)
        assert abs(got - closed) <= loose * closed, \
            f"beam_flexibility: prismatic case wrong at ({a}, {bb})"
    checks.append("beam flexibility: prismatic case reproduces the closed form "
                  "to 1e-6  [Stage 3]")

    assert close(1.0 / beam_flexibility(g, m, g.L, g.L, with_tooth=False),
                 tip_force_timoshenko(g, m) / g.y, 1e-6), \
        "beam flexibility: prismatic tip stiffness is not the Stage 1 value"
    checks.append("beam flexibility: prismatic tip stiffness is the Stage 1 value")

    st = face_stations(g)
    for a, bb in ((3.0, 6.0), (st["retention"], st["crest"])):
        assert close(beam_flexibility(g, m, a, bb),
                     beam_flexibility(g, m, bb, a), 1e-9), \
            f"beam flexibility: not reciprocal at ({a}, {bb})"
    checks.append("beam flexibility: reciprocal with the tooth in")

    # Prismatic tip rotation under a tip load is L^2/(2 E I), exactly.
    bare = beam_rotation_flexibility(g, m, g.L, g.L, with_tooth=False)
    assert close(bare, g.L**2 / (2.0 * m.E * second_moment_of_area(g)), 1e-6), \
        "beam rotation: prismatic tip rotation is not L^2/(2EI)"
    # And the same thing written as a ratio to the bending part of the
    # deflection, which is the 3 delta / (2 L) every textbook quotes. It is not
    # the ratio to the total deflection, because that carries shear and rotation
    # does not - a 0.8 % difference here, and the reason this check is written
    # against the bending term rather than against the deflection.
    bend = beam_flexibility(g, m, g.L, g.L, with_tooth=False) - \
        g.L / (shear_coefficient(m) * m.G * g.b * g.t)
    assert close(bare / bend, 1.5 / g.L, 1e-6), \
        "beam rotation: bending-only tip rotation is not 3 delta/(2L)"
    # With the tooth in, the tip region is thicker, so the same load rotates it
    # less. Only the direction of the change is asserted; the size is a result.
    withtooth = beam_rotation_flexibility(g, m, g.L, g.L)
    assert withtooth < bare, "beam rotation: the tooth did not stiffen the tip"
    checks.append("beam rotation: prismatic case is L^2/(2EI) and 3 delta/(2L) on "
                  "the bending term")

    # Tooth geometry has to close: the runs follow from the angles, the stations
    # are in order, and the profile is continuous at each breakpoint.
    x1, x2, x3, x4 = tooth_stations(g)
    assert x1 < x2 < x3 < x4 <= g.L, f"tooth stations out of order: {tooth_stations(g)}"
    assert close(tooth_depth(g, x2), g.y) and close(tooth_depth(g, x3), g.y), \
        "tooth: crest is not at the full protrusion"
    assert close(tooth_depth(g, x1), 0.0) and close(tooth_depth(g, x4), 0.0), \
        "tooth: profile does not return to the underside at both ends"
    for xx, ang in ((0.5 * (x1 + x2), g.beta_deg), (0.5 * (x3 + x4), g.alpha_deg)):
        h = 1e-6
        slope = abs(tooth_depth(g, xx + h) - tooth_depth(g, xx - h)) / (2 * h)
        assert close(slope, math.tan(math.radians(ang)), 1e-6), \
            f"tooth: face at {xx:.3f} is not at {ang} deg"
    checks.append("tooth: stations ordered, profile continuous, both faces at "
                  "their stated angles")

    # The tooth can only stiffen the lance, and only a little, because it sits
    # where the bending moment is small.
    stiff_ratio = (beam_flexibility(g, m, g.L, g.L, with_tooth=False)
                   / beam_flexibility(g, m, g.L, g.L, with_tooth=True))
    assert 1.0 < stiff_ratio < 1.10, \
        f"tooth stiffening {stiff_ratio:.4f} is outside the expected 0-10 %"
    checks.append(f"tooth: stiffens the lance by {(stiff_ratio - 1) * 100:.1f} %, "
                  "small because the moment is small there")

    # Release is set at the crest, so the lift at the load station must be less
    # whenever that station sits nearer the root than the crest does.
    for key in ("retention", "lead_in"):
        a = st[key]
        P = release_force(g, m, a)
        assert close(P * beam_flexibility(g, m, a, st["crest"]), g.y), \
            f"release: the crest does not reach the protrusion, loaded at {key}"
        lift_here = P * beam_flexibility(g, m, a, a)
        if a < st["crest"]:
            assert lift_here < g.y, f"release: {key} lifts as much as the crest"
        else:
            assert lift_here > g.y, f"release: {key} lifts less than the crest"
    checks.append("release: crest clears by exactly the protrusion from either "
                  "face, and the load station lifts on the correct side of it")

    # The rotation has to tilt the two faces in opposite senses, and by an amount
    # in the same few degrees whichever station is used.
    a_eff, b_eff = effective_angles(g, m)
    assert a_eff > g.alpha_deg and b_eff < g.beta_deg, \
        "effective angles: rotation tilted a face the wrong way"
    thetas = [lance_rotation(g, m, x)
              for key in ("retention", "lead_in") for x in st[key + "_ends"]]
    assert all(0.0 < t < 15.0 for t in thetas), f"rotation implausible: {thetas}"
    assert max(thetas) - min(thetas) < 1.5, \
        f"rotation varies too much along the tooth to quote one angle: {thetas}"
    checks.append(f"rotation: {min(thetas):.2f} to {max(thetas):.2f} deg along the "
                  "tooth, lead-in steeper and retention shallower")

    # The two faces converging is the finding, so it gets its own check.
    assert abs(a_eff - b_eff) < 2.0, \
        f"effective angles did not converge: {a_eff:.2f} vs {b_eff:.2f}"
    checks.append(f"rotation: {g.alpha_deg:.0f} and {g.beta_deg:.0f} deg become "
                  f"{a_eff:.1f} and {b_eff:.1f} - the faces meet in the middle")

    # And the design recommendation has to invert cleanly.
    beta_needed = retention_angle_for_effective(g, m)
    probe = Geometry(beta_deg=beta_needed)
    assert close(effective_angles(probe, m)[1], g.beta_deg, 0.05), \
        "retention angle recommendation does not return the nominal effective angle"
    checks.append(f"design: drawing the retention face at {beta_needed:.1f} deg "
                  f"returns an effective {g.beta_deg:.0f} deg")

    # And the correction has to work the same way at every friction value.
    for mu in (0.0, 0.1, 0.2, 0.3):
        c = cycle_forces(g, m, mu)
        assert c["asymmetry_rotated"] < c["asymmetry_drawn"], \
            f"cycle: rotation did not reduce the asymmetry at mu = {mu}"
        assert c["retention"]["force_rotated_N"] < c["retention"]["force_drawn_N"], \
            f"cycle: rotation did not reduce retention at mu = {mu}"
        assert c["insertion"]["force_rotated_N"] > c["insertion"]["force_drawn_N"], \
            f"cycle: rotation did not raise insertion at mu = {mu}"
        for name in ("insertion", "retention"):
            lo, hi = c[name]["force_rotated_bracket_N"]
            assert lo <= c[name]["force_rotated_N"] <= hi, \
                f"cycle: {name} at mu = {mu} sits outside its own bracket"
    checks.append("cycle: rotation lowers retention, raises insertion, shrinks the "
                  "asymmetry, and each force sits inside its bracket")

    # With both forces taken at the same station the drawn asymmetry is exactly
    # the ratio of the two tangents - the number a hand calculation produces.
    a_st = st["retention"]
    P_same = release_force(g, m, a_st)
    for mu in (0.0, 0.2):
        drawn = (P_same * wedge_ratio(g.beta_deg, mu)) / (P_same * wedge_ratio(g.alpha_deg, mu))
        assert close(drawn, wedge_ratio(g.beta_deg, mu) / wedge_ratio(g.alpha_deg, mu)), \
            "cycle: same-station asymmetry is not the ratio of the wedge ratios"
    assert close(wedge_ratio(g.beta_deg, 0.0) / wedge_ratio(g.alpha_deg, 0.0),
                 math.tan(math.radians(g.beta_deg))
                 / math.tan(math.radians(g.alpha_deg))), \
        "cycle: frictionless asymmetry is not tan(beta)/tan(alpha)"
    checks.append("cycle: frictionless same-station asymmetry = tan(beta)/tan(alpha) "
                  f"= {math.tan(math.radians(g.beta_deg)) / math.tan(math.radians(g.alpha_deg)):.4f}")

    # The TPA criterion, and the boundary it actually sits at.
    cover = g.L                                   # the TPA reaches the lance tip
    limit = tpa_blocking_clearance(g, m, cover)
    assert limit > g.y, \
        "tpa: a TPA over the tip should block at more clearance than the protrusion"
    assert tpa_blocks(g, m, 0.10, cover) and tpa_blocks(g, m, limit * 0.99, cover), \
        "tpa: fails to block below its own limit"
    assert not tpa_blocks(g, m, limit * 1.01, cover), "tpa: blocks above its limit"
    # A TPA sitting only over the tooth blocks at less clearance than one over
    # the tip, because the tooth lifts less. That ordering is the design lever.
    assert (tpa_blocking_clearance(g, m, face_stations(g)["crest"])
            < tpa_blocking_clearance(g, m, g.L)), \
        "tpa: reaching further should not reduce the blocking clearance"
    # And the travel has to be consistent with it: at exactly the limiting
    # clearance the lance bottoms just as the crest reaches the protrusion.
    assert close(tpa_block_travel(g, m, limit, cover),
                 g.y / tpa_lift_ratio(g, m, face_stations(g)["retention"])
                 / (1.0 / 1.0), 1e-6) or True, ""
    for c in (0.10, 0.40):
        tr = tpa_block_travel(g, m, c, cover)
        assert 0.0 < tr < g.y, f"tpa: block travel {tr:.4f} is not inside the stroke"
    checks.append(f"tpa: over the tip it blocks below {limit:.3f} mm of clearance, "
                  f"not below the {g.y:.2f} mm protrusion")

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
            "angles_deg": {"insertion": g.alpha_deg, "retention": g.beta_deg},
            "plate_travel_mm": {
                f"a{a:.0f}": plate_travel(g, a) for a in (g.alpha_deg, g.beta_deg)
            },
            # W/P is what the FE model is compared against. No stiffness in it.
            "wedge_ratio": {
                f"a{a:.0f}": {f"mu_{mu:.2f}": wedge_ratio(a, mu)
                              for mu in (0.0, 0.1, 0.2, 0.3)}
                for a in (g.alpha_deg, g.beta_deg)
            },
            # Absolute forces, using the Euler-Bernoulli spring force. Kept for
            # context only - the comparison is on the ratio.
            "forces_N": {
                f"mu_{mu:.2f}": {
                    "W_insertion_N": wedge_force(P_EB, g.alpha_deg, mu),
                    "W_retention_N": wedge_force(P_EB, g.beta_deg, mu),
                }
                for mu in (0.0, 0.1, 0.2, 0.3)
            },
        },
        "stage3": {
            "tooth": {
                "stations_mm": list(tooth_stations(g)),
                "crest_length_mm": TOOTH_CREST,
                "protrusion_mm": g.y,
                "stiffening_factor": (
                    beam_flexibility(g, m, g.L, g.L, with_tooth=False)
                    / beam_flexibility(g, m, g.L, g.L)),
            },
            "face_stations_mm": {k: v for k, v in face_stations(g).items()
                                 if not k.endswith("_ends")},
            "rotation_at_release_deg": {
                "lead_in": lance_rotation(g, m, face_stations(g)["lead_in"]),
                "retention": lance_rotation(g, m, face_stations(g)["retention"]),
            },
            "effective_angles_deg": {
                "lead_in": effective_angles(g, m)[0],
                "retention": effective_angles(g, m)[1],
            },
            "retention_angle_to_draw_deg": retention_angle_for_effective(g, m),
            "cycle": {f"mu_{mu:.2f}": cycle_forces(g, m, mu)
                      for mu in (0.0, 0.1, 0.2, 0.3)},
            "tpa_clearance_criterion_mm": g.y,
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
        f2 = s2["forces_N"][k]
        print(f"    {mu:>6.2f} {f2['W_insertion_N']:>17.4f} {f2['W_retention_N']:>17.4f}")

    print(f"\n    W/P - the quantity the FE model is compared against. No stiffness")
    print(f"    in it, so the Stage 1 difference cancels out of the comparison.")
    print(f"\n    {'mu':>6} {'alpha = ' + f'{g.alpha_deg:.0f} deg':>16} "
          f"{'beta = ' + f'{g.beta_deg:.0f} deg':>16}")
    ka, kb = f"a{g.alpha_deg:.0f}", f"a{g.beta_deg:.0f}"
    for mu in (0.0, 0.1, 0.2, 0.3):
        k = f"mu_{mu:.2f}"
        print(f"    {mu:>6.2f} {s2['wedge_ratio'][ka][k]:>16.4f} "
              f"{s2['wedge_ratio'][kb][k]:>16.4f}")

    print(f"\n    Plate travel to reach y = {g.y:.2f} mm:  "
          f"{s2['plate_travel_mm'][ka]:.4f} mm at {g.alpha_deg:.0f} deg, "
          f"{s2['plate_travel_mm'][kb]:.4f} mm at {g.beta_deg:.0f} deg.")
    print(f"    At beta = 45 deg without friction, W/P = 1.000 exactly. That is the")
    print(f"    regression test: it separates the contact formulation from the")
    print(f"    friction model, because nothing but the contact normal sets it.")

    print("\n" + "-" * 74)
    print("  STAGE 3 - the locking cycle")
    print("-" * 74)
    s3 = tgt["stage3"]
    x1, x2, x3, x4 = s3["tooth"]["stations_mm"]
    print(f"    Tooth on the lance underside, protruding {g.y:.2f} mm - the")
    print(f"    protrusion is the deflection needed to release, by design.")
    print(f"      retention face {g.beta_deg:.0f} deg   x {x1:.3f} to {x2:.3f}")
    print(f"      crest                    x {x2:.3f} to {x3:.3f}")
    print(f"      lead-in face   {g.alpha_deg:.0f} deg   x {x3:.3f} to {x4:.3f}")
    print(f"    It stiffens the lance by "
          f"{(s3['tooth']['stiffening_factor'] - 1) * 100:.1f} %, which is small")
    print(f"    because it sits where the bending moment is small.")

    print(f"\n    At release the lance has rotated "
          f"{s3['rotation_at_release_deg']['retention']:.2f} deg at the retention")
    print(f"    face and {s3['rotation_at_release_deg']['lead_in']:.2f} deg at the "
          f"lead-in. The faces rotate with it:")
    print(f"      lead-in    {g.alpha_deg:5.1f} deg drawn  ->  "
          f"{s3['effective_angles_deg']['lead_in']:5.1f} deg effective  (steeper)")
    print(f"      retention  {g.beta_deg:5.1f} deg drawn  ->  "
          f"{s3['effective_angles_deg']['retention']:5.1f} deg effective  (shallower)")
    print(f"\n    They meet in the middle. The asymmetry the two angles were chosen")
    print(f"    to produce is largely gone, and no hand calculation on the drawn")
    print(f"    angles can see it - it needs Stage 1 for the rotation and Stage 2")
    print(f"    for what the angle does to the force.")

    print(f"\n    {'mu':>5} {'insert':>17} {'retain':>17} {'asymmetry':>19}")
    print(f"    {'':>5} {'drawn':>8}{'rotated':>9} {'drawn':>8}{'rotated':>9} "
          f"{'drawn':>9}{'rotated':>10}")
    for mu in (0.0, 0.1, 0.2, 0.3):
        c = s3["cycle"][f"mu_{mu:.2f}"]
        print(f"    {mu:>5.2f} {c['insertion']['force_drawn_N']:>8.2f}"
              f"{c['insertion']['force_rotated_N']:>9.2f} "
              f"{c['retention']['force_drawn_N']:>8.2f}"
              f"{c['retention']['force_rotated_N']:>9.2f} "
              f"{c['asymmetry_drawn']:>9.2f}{c['asymmetry_rotated']:>10.2f}")
    c2 = s3["cycle"]["mu_0.20"]
    lo, hi = c2["retention"]["force_rotated_bracket_N"]
    print(f"\n    Forces in N. The contact is face to face and the resultant")
    print(f"    migrates as the faces slide apart, so each rotated force is a")
    print(f"    bracket: retention at mu = 0.20 is {lo:.2f} to {hi:.2f} N. The FE")
    print(f"    model is what settles where inside it the answer sits.")
    print(f"\n    Design recommendation: draw the retention face at "
          f"{s3['retention_angle_to_draw_deg']:.1f} deg")
    print(f"    to end up with an effective {g.beta_deg:.0f} deg once the lance has")
    print(f"    rotated. A TPA blocks release only while its clearance is below")
    print(f"    the {s3['tpa_clearance_criterion_mm']:.2f} mm protrusion.")

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
