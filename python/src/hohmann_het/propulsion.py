# Author: A Taylor | Purpose: Hall Effect Thruster propulsion models | Ref: Vallado/Curtis
"""High-fidelity Hall Effect Thruster (HET) propulsion module.

Models anode efficiency, beam voltage, exhaust velocity, specific impulse,
mass flow rate, and propellant budget.

References
----------
- Goebel & Katz, "Fundamentals of Electric Propulsion," JPL/Wiley, 2008.
- Vallado, "Fundamentals of Astrodynamics and Applications," 4th ed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import astropy.units as u

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
#: Standard gravity [m/s^2]
G0: u.Quantity = 9.80665 * u.m / u.s**2

#: Elementary charge [C]
E_CHARGE: float = 1.602176634e-19

#: Xenon atom mass [kg]  (131.293 u, 1 u = 1.66053906660e-27 kg)
M_XE: float = 2.180174e-25


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HETState:
    """Hall Effect Thruster operating-point state.

    Attributes
    ----------
    isp : u.Quantity
        Specific impulse [s].
    exhaust_velocity : u.Quantity
        Effective exhaust velocity [m/s].
    thrust : u.Quantity
        Thrust force [N].
    mass_flow : u.Quantity
        Propellant mass flow rate [kg/s].
    anode_efficiency : float
        Anode efficiency (dimensionless).
    discharge_power : u.Quantity
        Discharge power [W].
    """

    isp: u.Quantity
    exhaust_velocity: u.Quantity
    thrust: u.Quantity
    mass_flow: u.Quantity
    anode_efficiency: float
    discharge_power: u.Quantity


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------
def compute_het_state(
    discharge_voltage: "u.Quantity | float",
    discharge_power: "u.Quantity | float",
    anode_efficiency: float,
) -> HETState:
    r"""Compute Hall Effect Thruster operating point from discharge parameters.

    **Beam-limited exhaust velocity** (Goebel & Katz, Eq. 7-13):

    .. math::
        v_e = \sqrt{\frac{2 \eta_a e V_d}{m_{\rm Xe}}}

    **Specific impulse:**

    .. math::
        I_{sp} = \frac{v_e}{g_0}

    **Thrust** (power balance):

    .. math::
        F = \frac{2 \eta_a P_d}{v_e}

    **Mass flow rate:**

    .. math::
        \dot{m} = \frac{F}{v_e}

    Parameters
    ----------
    discharge_voltage : astropy.units.Quantity or float
        Anode discharge voltage [V]. If plain ``float``, assumed V.
    discharge_power : astropy.units.Quantity or float
        Input discharge power [W]. If plain ``float``, assumed W.
    anode_efficiency : float
        Anode efficiency :math:`\eta_a \in (0, 1)`.

    Returns
    -------
    HETState
        Complete HET operating point with units.

    Raises
    ------
    ValueError
        If voltage, power, or efficiency are out of range.
    """
    if not isinstance(discharge_voltage, u.Quantity):
        discharge_voltage = float(discharge_voltage) * u.V
    if not isinstance(discharge_power, u.Quantity):
        discharge_power = float(discharge_power) * u.W

    V_d = discharge_voltage.to(u.V).value
    P_d = discharge_power.to(u.W).value

    if V_d <= 0.0 or P_d <= 0.0:
        raise ValueError("Discharge voltage and power must be positive.")
    if not (0.0 < anode_efficiency < 1.0):
        raise ValueError("Anode efficiency must be in (0, 1).")

    ve   = math.sqrt(2.0 * anode_efficiency * E_CHARGE * V_d / M_XE)  # m/s
    isp  = ve / G0.value                                                # s
    F    = 2.0 * anode_efficiency * P_d / ve                           # N
    mdot = F / ve                                                       # kg/s

    return HETState(
        isp=isp * u.s,
        exhaust_velocity=ve * (u.m / u.s),
        thrust=F * u.N,
        mass_flow=mdot * (u.kg / u.s),
        anode_efficiency=anode_efficiency,
        discharge_power=P_d * u.W,
    )


def propellant_mass(
    m_initial: "u.Quantity | float",
    delta_v: "u.Quantity | float",
    isp: "u.Quantity | float",
) -> u.Quantity:
    r"""Propellant mass from the Tsiolkovsky rocket equation.

    .. math::
        m_p = m_0 \left(1 - \exp\!\left(-\frac{\Delta v}{g_0 \, I_{sp}}\right)\right)

    Parameters
    ----------
    m_initial : astropy.units.Quantity or float
        Initial (wet) spacecraft mass [kg]. If plain ``float``, assumed kg.
    delta_v : astropy.units.Quantity or float
        Required delta-V. If plain ``float``, assumed km/s.
    isp : astropy.units.Quantity or float
        Specific impulse [s]. If plain ``float``, assumed s.

    Returns
    -------
    u.Quantity
        Propellant mass [kg].
    """
    if not isinstance(m_initial, u.Quantity):
        m_initial = float(m_initial) * u.kg
    if not isinstance(delta_v, u.Quantity):
        delta_v = float(delta_v) * (u.km / u.s)
    if not isinstance(isp, u.Quantity):
        isp = float(isp) * u.s

    m0   = m_initial.to(u.kg).value
    dv   = delta_v.to(u.m / u.s).value   # convert km/s -> m/s
    isp_s = isp.to(u.s).value

    ve = isp_s * G0.value                 # m/s
    mp = m0 * (1.0 - math.exp(-dv / ve))
    return mp * u.kg


def burn_time(
    thrust: "u.Quantity | float",
    mass_flow: "u.Quantity | float",
    delta_v: "u.Quantity | float",
    m_initial: "u.Quantity | float",
) -> u.Quantity:
    r"""Constant-thrust burn time for a given delta-V.

    Derived from integrating the rocket equation at constant mass flow:

    .. math::
        T_b = \frac{m_0}{\dot{m}}
              \left(1 - \exp\!\left(-\frac{\Delta v}{v_e}\right)\right)

    where :math:`v_e = F / \dot{m}`.

    Parameters
    ----------
    thrust : astropy.units.Quantity or float
        Thruster force [N]. If plain ``float``, assumed N.
    mass_flow : astropy.units.Quantity or float
        Mass flow rate [kg/s]. If plain ``float``, assumed kg/s.
    delta_v : astropy.units.Quantity or float
        Required delta-V. If plain ``float``, assumed km/s.
    m_initial : astropy.units.Quantity or float
        Initial spacecraft mass [kg]. If plain ``float``, assumed kg.

    Returns
    -------
    u.Quantity
        Burn time [s].
    """
    if not isinstance(thrust, u.Quantity):
        thrust = float(thrust) * u.N
    if not isinstance(mass_flow, u.Quantity):
        mass_flow = float(mass_flow) * (u.kg / u.s)
    if not isinstance(delta_v, u.Quantity):
        delta_v = float(delta_v) * (u.km / u.s)
    if not isinstance(m_initial, u.Quantity):
        m_initial = float(m_initial) * u.kg

    F    = thrust.to(u.N).value
    mdot = mass_flow.to(u.kg / u.s).value
    dv   = delta_v.to(u.m / u.s).value
    m0   = m_initial.to(u.kg).value

    if F <= 0.0 or mdot <= 0.0 or m0 <= 0.0 or dv < 0.0:
        raise ValueError("Thrust, mass flow, and initial mass must be positive; dv >= 0.")

    ve = F / mdot                                    # m/s
    tb = (m0 / mdot) * (1.0 - math.exp(-dv / ve))
    return tb * u.s
