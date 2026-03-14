# Author: A Taylor | Purpose: Mission optimization solvers for HET transfers | Ref: Vallado/Curtis
"""Mission optimization module.

Provides solvers for minimizing propellant mass or burn time for
low-thrust Hall Effect Thruster orbital transfer missions.

The core algorithm is a pure-Python golden-section search, ensuring
identical numerical behavior across the tri-language implementation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import astropy.units as u

from .propulsion import G0

# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class OptimizationResult:
    """Result of a specific-impulse optimization.

    Attributes
    ----------
    optimal_isp : u.Quantity
        Optimal specific impulse [s].
    propellant_mass : u.Quantity
        Propellant mass at the optimum [kg].
    burn_time : u.Quantity
        Burn time at the optimum [s].
    objective_value : float
        Dimensionless composite objective at the optimum.
    """

    optimal_isp: u.Quantity
    propellant_mass: u.Quantity
    burn_time: u.Quantity
    objective_value: float


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------
def golden_section_minimize(
    func: Callable[[float], float],
    a: float,
    b: float,
    tol: float = 1e-9,
) -> float:
    r"""Golden-section search for the minimum of a unimodal scalar function.

    Converges to within ``tol`` of the true minimum location.

    The interval is reduced at each step by the golden ratio:

    .. math::
        \varphi = \frac{1 + \sqrt{5}}{2} \approx 1.618

    Parameters
    ----------
    func : callable
        Scalar function :math:`f : \mathbb{R} \to \mathbb{R}` to minimize.
    a : float
        Left bound of the search interval.
    b : float
        Right bound of the search interval.
    tol : float
        Convergence tolerance on the interval width.

    Returns
    -------
    float
        Approximate location of the minimum.
    """
    gr = (math.sqrt(5.0) + 1.0) / 2.0  # golden ratio
    c = b - (b - a) / gr
    d = a + (b - a) / gr

    while abs(b - a) > tol:
        if func(c) < func(d):
            b = d
        else:
            a = c
        c = b - (b - a) / gr
        d = a + (b - a) / gr

    return (a + b) / 2.0


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def optimize_isp(
    m_initial: "u.Quantity | float",
    delta_v: "u.Quantity | float",
    discharge_power: "u.Quantity | float",
    anode_efficiency: float,
    isp_min: "u.Quantity | float" = 500.0,
    isp_max: "u.Quantity | float" = 5000.0,
    lambda_weight: float = 1e-4,
) -> OptimizationResult:
    r"""Optimize specific impulse to minimize a composite propellant-time objective.

    **Objective function** (Pareto trade-off between propellant mass and burn time):

    .. math::
        J(I_{sp}) = \frac{m_p(I_{sp})}{m_0}
                  + \lambda \frac{T_b(I_{sp})}{T_{\rm ref}}

    **Propellant mass** (Tsiolkovsky):

    .. math::
        m_p = m_0 \left(1 - e^{-\Delta v / (g_0 I_{sp})}\right)

    **Power-limited burn time:**

    .. math::
        T_b = \frac{m_0 \,\Delta v \, v_e}{2 \eta_a P_d}, \quad v_e = g_0 I_{sp}

    Minimized via :func:`golden_section_minimize` over
    :math:`I_{sp} \in [I_{sp,\min},\ I_{sp,\max}]`.

    Parameters
    ----------
    m_initial : astropy.units.Quantity or float
        Spacecraft initial (wet) mass [kg].
    delta_v : astropy.units.Quantity or float
        Required delta-V. If plain ``float``, assumed km/s.
    discharge_power : astropy.units.Quantity or float
        Available thruster discharge power [W].
    anode_efficiency : float
        HET anode efficiency :math:`\eta_a \in (0, 1)`.
    isp_min : astropy.units.Quantity or float
        Lower bound of Isp search [s]. Default 500 s.
    isp_max : astropy.units.Quantity or float
        Upper bound of Isp search [s]. Default 5000 s.
    lambda_weight : float
        Dimensionless weighting factor for the burn-time term.

    Returns
    -------
    OptimizationResult
        Optimal Isp and corresponding mission parameters.
    """
    if not isinstance(m_initial, u.Quantity):
        m_initial = float(m_initial) * u.kg
    if not isinstance(delta_v, u.Quantity):
        delta_v = float(delta_v) * (u.km / u.s)
    if not isinstance(discharge_power, u.Quantity):
        discharge_power = float(discharge_power) * u.W
    if not isinstance(isp_min, u.Quantity):
        isp_min = float(isp_min) * u.s
    if not isinstance(isp_max, u.Quantity):
        isp_max = float(isp_max) * u.s

    m0     = m_initial.to(u.kg).value
    dv_ms  = delta_v.to(u.m / u.s).value
    P      = discharge_power.to(u.W).value
    isp_lo = isp_min.to(u.s).value
    isp_hi = isp_max.to(u.s).value
    eta_a  = anode_efficiency
    g0     = G0.value

    if m0 <= 0.0 or dv_ms < 0.0 or P <= 0.0:
        raise ValueError("Mass, delta-V, and power must be non-negative/positive.")
    if not (0.0 < eta_a < 1.0):
        raise ValueError("Anode efficiency must be in (0, 1).")
    if isp_lo <= 0.0 or isp_hi <= isp_lo:
        raise ValueError("Isp bounds must satisfy 0 < isp_min < isp_max.")

    # Normalisation reference: burn time at minimum Isp
    ve_ref = isp_lo * g0
    T_ref  = m0 * dv_ms * ve_ref / (2.0 * eta_a * P)

    def objective(isp_s: float) -> float:
        ve      = isp_s * g0
        mp_frac = 1.0 - math.exp(-dv_ms / ve)
        t_burn  = m0 * dv_ms * ve / (2.0 * eta_a * P)
        return mp_frac + lambda_weight * (t_burn / T_ref)

    isp_opt = golden_section_minimize(objective, isp_lo, isp_hi)
    ve_opt  = isp_opt * g0
    mp_opt  = m0 * (1.0 - math.exp(-dv_ms / ve_opt))
    tb_opt  = m0 * dv_ms * ve_opt / (2.0 * eta_a * P)
    J_opt   = objective(isp_opt)

    return OptimizationResult(
        optimal_isp=isp_opt * u.s,
        propellant_mass=mp_opt * u.kg,
        burn_time=tb_opt * u.s,
        objective_value=J_opt,
    )


def min_propellant_transfer(
    h1: "u.Quantity | float",
    h2: "u.Quantity | float",
    m_initial: "u.Quantity | float",
    discharge_power: "u.Quantity | float",
    anode_efficiency: float,
    isp_min: "u.Quantity | float" = 500.0,
    isp_max: "u.Quantity | float" = 5000.0,
) -> OptimizationResult:
    """Compute Hohmann delta-V then optimize Isp for minimum propellant.

    Convenience wrapper combining :func:`~hohmann_het.dynamics.compute_hohmann`
    and :func:`optimize_isp`.

    Parameters
    ----------
    h1 : astropy.units.Quantity or float
        Initial orbit altitude [km].
    h2 : astropy.units.Quantity or float
        Final orbit altitude [km].
    m_initial : astropy.units.Quantity or float
        Spacecraft wet mass [kg].
    discharge_power : astropy.units.Quantity or float
        Available HET discharge power [W].
    anode_efficiency : float
        Anode efficiency in (0, 1).
    isp_min : astropy.units.Quantity or float
        Lower Isp search bound [s].
    isp_max : astropy.units.Quantity or float
        Upper Isp search bound [s].

    Returns
    -------
    OptimizationResult
    """
    from .dynamics import compute_hohmann

    transfer = compute_hohmann(h1, h2)
    return optimize_isp(
        m_initial=m_initial,
        delta_v=transfer.total_dv,
        discharge_power=discharge_power,
        anode_efficiency=anode_efficiency,
        isp_min=isp_min,
        isp_max=isp_max,
    )
