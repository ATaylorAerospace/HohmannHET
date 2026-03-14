// Author: A Taylor | Purpose: Hall Effect Thruster propulsion models | Ref: Vallado/Curtis
#pragma once

/**
 * @file propulsion.hpp
 * @brief High-fidelity Hall Effect Thruster (HET) propulsion models.
 *
 * Models anode efficiency, beam voltage, exhaust velocity, specific impulse,
 * mass flow rate, and propellant budget.
 *
 * Reference:
 *   Goebel & Katz, "Fundamentals of Electric Propulsion," JPL/Wiley, 2008.
 *
 * @author A Taylor
 * @ref Vallado/Curtis
 */

#include <cmath>
#include <stdexcept>

#include "dynamics.hpp"   // for constants::G0

namespace hohmann_het {

// ---------------------------------------------------------------------------
// Propulsion-specific physical constants
// ---------------------------------------------------------------------------
namespace prop_constants {

/// Elementary charge [C]
inline constexpr double E_CHARGE = 1.602176634e-19;

/// Xenon atom mass [kg]  (131.293 u, 1 u = 1.66053906660e-27 kg)
inline constexpr double M_XE = 2.180174e-25;

} // namespace prop_constants


// ---------------------------------------------------------------------------
// Result struct
// ---------------------------------------------------------------------------
/**
 * @brief HET operating-point state variables.
 */
struct HETState {
    double isp;               ///< Specific impulse [s]
    double exhaust_velocity;  ///< Effective exhaust velocity [m/s]
    double thrust;            ///< Thrust force [N]
    double mass_flow;         ///< Propellant mass flow rate [kg/s]
    double anode_efficiency;  ///< Anode efficiency [-]
    double discharge_power;   ///< Discharge power [W]
};


// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------
/**
 * @brief Compute HET operating point from discharge parameters.
 *
 * @par Mathematical formulation
 *
 * Beam-limited exhaust velocity (Goebel & Katz, Eq. 7-13):
 * @f[ v_e = \sqrt{\frac{2 \eta_a e V_d}{m_{\rm Xe}}} @f]
 *
 * Specific impulse:
 * @f[ I_{sp} = \frac{v_e}{g_0} @f]
 *
 * Thrust (power balance):
 * @f[ F = \frac{2 \eta_a P_d}{v_e} @f]
 *
 * Mass flow rate:
 * @f[ \dot{m} = \frac{F}{v_e} @f]
 *
 * @param discharge_voltage  Anode discharge voltage [V]. Must be positive.
 * @param discharge_power    Input discharge power [W]. Must be positive.
 * @param anode_efficiency   Anode efficiency @f$\eta_a \in (0, 1)@f$.
 * @return                   HETState with complete operating point.
 * @throws std::invalid_argument on invalid inputs.
 */
[[nodiscard]] inline HETState compute_het_state(
    double discharge_voltage,
    double discharge_power,
    double anode_efficiency)
{
    using namespace prop_constants;

    if (discharge_voltage <= 0.0 || discharge_power <= 0.0)
        throw std::invalid_argument("Discharge voltage and power must be positive.");
    if (anode_efficiency <= 0.0 || anode_efficiency >= 1.0)
        throw std::invalid_argument("Anode efficiency must be in (0, 1).");

    const double ve    = std::sqrt(2.0 * anode_efficiency * E_CHARGE
                                   * discharge_voltage / M_XE);
    const double isp   = ve / constants::G0;
    const double F     = 2.0 * anode_efficiency * discharge_power / ve;
    const double mdot  = F / ve;

    return HETState{isp, ve, F, mdot, anode_efficiency, discharge_power};
}


/**
 * @brief Propellant mass from the Tsiolkovsky rocket equation.
 *
 * @f[
 *   m_p = m_0 \left(1 - \exp\!\left(-\frac{\Delta v}{g_0 \, I_{sp}}\right)\right)
 * @f]
 *
 * @param m_initial  Initial (wet) spacecraft mass [kg].
 * @param delta_v    Required delta-V [m/s]. Non-negative.
 * @param isp        Specific impulse [s]. Must be positive.
 * @return           Propellant mass [kg].
 * @throws std::invalid_argument on invalid inputs.
 */
[[nodiscard]] inline double propellant_mass(
    double m_initial,
    double delta_v,
    double isp)
{
    if (m_initial <= 0.0 || isp <= 0.0)
        throw std::invalid_argument("m_initial and isp must be positive.");
    if (delta_v < 0.0)
        throw std::invalid_argument("delta_v must be non-negative.");

    const double ve = isp * constants::G0;
    return m_initial * (1.0 - std::exp(-delta_v / ve));
}


/**
 * @brief Constant-thrust burn time for a given delta-V.
 *
 * @f[
 *   T_b = \frac{m_0}{\dot{m}}
 *         \left(1 - \exp\!\left(-\frac{\Delta v}{v_e}\right)\right), \quad
 *   v_e = \frac{F}{\dot{m}}
 * @f]
 *
 * @param thrust     Thruster force [N]. Must be positive.
 * @param mass_flow  Mass flow rate [kg/s]. Must be positive.
 * @param delta_v    Required delta-V [m/s]. Non-negative.
 * @param m_initial  Initial spacecraft mass [kg]. Must be positive.
 * @return           Burn time [s].
 * @throws std::invalid_argument on invalid inputs.
 */
[[nodiscard]] inline double burn_time(
    double thrust,
    double mass_flow,
    double delta_v,
    double m_initial)
{
    if (thrust <= 0.0 || mass_flow <= 0.0 || m_initial <= 0.0)
        throw std::invalid_argument("Thrust, mass_flow, and m_initial must be positive.");
    if (delta_v < 0.0)
        throw std::invalid_argument("delta_v must be non-negative.");

    const double ve = thrust / mass_flow;
    return (m_initial / mass_flow) * (1.0 - std::exp(-delta_v / ve));
}

} // namespace hohmann_het
