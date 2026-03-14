# Author: A Taylor | Purpose: Keplerian dynamics and Hohmann transfer | Ref: Vallado/Curtis
"""Keplerian orbital dynamics module.

Implements two-impulse Hohmann transfer orbit computations following:
  - Vallado, "Fundamentals of Astrodynamics and Applications," 4th ed.
  - Curtis, "Orbital Mechanics for Engineering Students," 3rd ed.

All physical quantities are represented as astropy.units.Quantity objects.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import astropy.units as u

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
#: Earth gravitational parameter [km^3 s^-2]
MU: u.Quantity = 398600.4418 * u.km**3 / u.s**2

#: Earth equatorial radius [km]
R_EARTH: u.Quantity = 6378.137 * u.km


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HohmannTransfer:
    """Results of a two-impulse Hohmann transfer computation.

    All attributes are :class:`astropy.units.Quantity` objects carrying
    explicit SI/km units.

    Attributes
    ----------
    r1 : u.Quantity
        Initial circular orbit radius [km].
    r2 : u.Quantity
        Final circular orbit radius [km].
    a_transfer : u.Quantity
        Transfer ellipse semi-major axis [km].
    dv_departure : u.Quantity
        First impulsive delta-V (perigee burn) [km/s].
    dv_arrival : u.Quantity
        Second impulsive delta-V (apogee burn) [km/s].
    tof : u.Quantity
        Transfer time-of-flight [s].
    """

    r1: u.Quantity
    r2: u.Quantity
    a_transfer: u.Quantity
    dv_departure: u.Quantity
    dv_arrival: u.Quantity
    tof: u.Quantity

    @property
    def total_dv(self) -> u.Quantity:
        r"""Total delta-V for the transfer.

        .. math::
            \Delta v_{\rm total} = \left|\Delta v_1\right| + \left|\Delta v_2\right|

        Returns
        -------
        u.Quantity
            Total delta-V [km/s].
        """
        return abs(self.dv_departure) + abs(self.dv_arrival)

    @property
    def tof_hours(self) -> u.Quantity:
        """Transfer time in hours."""
        return self.tof.to(u.hour)


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------
def compute_hohmann(
    h1: "u.Quantity | float",
    h2: "u.Quantity | float",
) -> HohmannTransfer:
    r"""Compute a two-impulse Hohmann transfer between circular orbits.

    Follows Vallado Section 6.3 (impulsive Hohmann transfer).

    **Orbital radii:**

    .. math::
        r_i = R_\oplus + h_i

    **Transfer semi-major axis:**

    .. math::
        a_t = \frac{r_1 + r_2}{2}

    **Circular orbit velocities:**

    .. math::
        v_{c,i} = \sqrt{\frac{\mu}{r_i}}

    **Transfer orbit velocities at periapsis / apoapsis:**

    .. math::
        v_{t,1} = \sqrt{\mu \left(\frac{2}{r_1} - \frac{1}{a_t}\right)}, \quad
        v_{t,2} = \sqrt{\mu \left(\frac{2}{r_2} - \frac{1}{a_t}\right)}

    **Impulsive burns:**

    .. math::
        \Delta v_1 = v_{t,1} - v_{c,1}, \quad
        \Delta v_2 = v_{c,2} - v_{t,2}

    **Transfer time-of-flight (half-period of transfer ellipse):**

    .. math::
        T = \pi \sqrt{\frac{a_t^3}{\mu}}

    Parameters
    ----------
    h1 : astropy.units.Quantity or float
        Initial orbit altitude above Earth's surface.
        If plain ``float``, interpreted as km.
    h2 : astropy.units.Quantity or float
        Final orbit altitude above Earth's surface.
        If plain ``float``, interpreted as km.

    Returns
    -------
    HohmannTransfer
        Dataclass with all transfer parameters.

    Raises
    ------
    ValueError
        If either altitude is non-positive.

    Examples
    --------
    >>> import astropy.units as u
    >>> from hohmann_het.dynamics import compute_hohmann
    >>> r = compute_hohmann(400.0 * u.km, 35786.0 * u.km)
    >>> print(f"total dv = {r.total_dv:.4f}")
    total dv = 3.8526 km / s
    """
    if not isinstance(h1, u.Quantity):
        h1 = float(h1) * u.km
    if not isinstance(h2, u.Quantity):
        h2 = float(h2) * u.km

    h1 = h1.to(u.km)
    h2 = h2.to(u.km)

    if h1.value <= 0.0 or h2.value <= 0.0:
        raise ValueError("Orbit altitudes must be positive.")

    r1 = R_EARTH + h1
    r2 = R_EARTH + h2
    a_t = (r1 + r2) * 0.5

    # Work in raw floats for speed; attach units to results
    mu = MU.value        # km^3/s^2
    r1v = r1.value       # km
    r2v = r2.value       # km
    av  = a_t.value      # km

    mu_over_a = mu / av

    vc1 = math.sqrt(mu / r1v)
    vc2 = math.sqrt(mu / r2v)
    vt1 = math.sqrt(2.0 * mu / r1v - mu_over_a)
    vt2 = math.sqrt(2.0 * mu / r2v - mu_over_a)

    dv1 = (vt1 - vc1) * (u.km / u.s)
    dv2 = (vc2 - vt2) * (u.km / u.s)
    tof = math.pi * av * math.sqrt(av / mu) * u.s

    return HohmannTransfer(
        r1=r1,
        r2=r2,
        a_transfer=a_t,
        dv_departure=dv1,
        dv_arrival=dv2,
        tof=tof,
    )


def circular_velocity(r: "u.Quantity | float") -> u.Quantity:
    r"""Circular orbit speed at radius *r*.

    .. math::
        v_c = \sqrt{\frac{\mu}{r}}

    Parameters
    ----------
    r : astropy.units.Quantity or float
        Orbital radius. If plain ``float``, interpreted as km.

    Returns
    -------
    u.Quantity
        Circular velocity [km/s].
    """
    if not isinstance(r, u.Quantity):
        r = float(r) * u.km
    r = r.to(u.km)
    if r.value <= 0.0:
        raise ValueError("Orbital radius must be positive.")
    return math.sqrt(MU.value / r.value) * (u.km / u.s)


def orbital_period(a: "u.Quantity | float") -> u.Quantity:
    r"""Keplerian orbital period.

    .. math::
        T = 2\pi \sqrt{\frac{a^3}{\mu}}

    Parameters
    ----------
    a : astropy.units.Quantity or float
        Semi-major axis. If plain ``float``, interpreted as km.

    Returns
    -------
    u.Quantity
        Orbital period [s].
    """
    if not isinstance(a, u.Quantity):
        a = float(a) * u.km
    a = a.to(u.km)
    if a.value <= 0.0:
        raise ValueError("Semi-major axis must be positive.")
    av = a.value
    return 2.0 * math.pi * av * math.sqrt(av / MU.value) * u.s
