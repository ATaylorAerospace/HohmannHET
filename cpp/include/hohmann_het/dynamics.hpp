// Author: A Taylor | Purpose: Keplerian dynamics and Hohmann transfer | Ref: Vallado/Curtis
#pragma once

/**
 * @file dynamics.hpp
 * @brief Keplerian orbital dynamics and two-impulse Hohmann transfer.
 *
 * Implements the Hohmann transfer following:
 *   - Vallado, "Fundamentals of Astrodynamics and Applications," 4th ed.
 *   - Curtis, "Orbital Mechanics for Engineering Students," 3rd ed.
 *
 * All lengths in km, velocities in km/s, time in seconds unless noted.
 *
 * @author A Taylor
 * @ref Vallado/Curtis
 */

#include <cmath>
#include <stdexcept>

namespace hohmann_het {

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------
namespace constants {

/// Earth gravitational parameter [km^3/s^2]
inline constexpr double MU = 398600.4418;

/// Earth equatorial radius [km]
inline constexpr double R_EARTH = 6378.137;

/// Standard gravity [m/s^2]
inline constexpr double G0 = 9.80665;

/// Mathematical constant pi
inline constexpr double PI = 3.14159265358979323846;

} // namespace constants


// ---------------------------------------------------------------------------
// Result struct
// ---------------------------------------------------------------------------
/**
 * @brief Results of a two-impulse Hohmann transfer computation.
 *
 * All lengths in km, velocities in km/s, time in seconds.
 */
struct HohmannTransfer {
    double r1;            ///< Initial orbit radius [km]
    double r2;            ///< Final   orbit radius [km]
    double a_transfer;    ///< Transfer ellipse semi-major axis [km]
    double dv_departure;  ///< First  impulsive delta-V [km/s]
    double dv_arrival;    ///< Second impulsive delta-V [km/s]
    double tof;           ///< Transfer time-of-flight [s]

    /// Total delta-V (absolute sum of both burns) [km/s]
    [[nodiscard]] double total_dv() const noexcept {
        return std::abs(dv_departure) + std::abs(dv_arrival);
    }

    /// Transfer time-of-flight in hours
    [[nodiscard]] constexpr double tof_hours() const noexcept {
        return tof / 3600.0;
    }
};


// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------
/**
 * @brief Compute a two-impulse Hohmann transfer between circular orbits.
 *
 * @par Mathematical formulation
 *
 * Orbital radii:
 * @f[ r_i = R_\oplus + h_i @f]
 *
 * Transfer semi-major axis:
 * @f[ a_t = \frac{r_1 + r_2}{2} @f]
 *
 * Circular orbit velocities:
 * @f[ v_{c,i} = \sqrt{\frac{\mu}{r_i}} @f]
 *
 * Transfer orbit velocities at periapsis / apoapsis:
 * @f[ v_{t,1} = \sqrt{\mu \left(\frac{2}{r_1} - \frac{1}{a_t}\right)}, \quad
 *     v_{t,2} = \sqrt{\mu \left(\frac{2}{r_2} - \frac{1}{a_t}\right)} @f]
 *
 * Impulsive burns:
 * @f[ \Delta v_1 = v_{t,1} - v_{c,1}, \quad
 *     \Delta v_2 = v_{c,2} - v_{t,2} @f]
 *
 * Transfer time-of-flight (half-period of transfer ellipse):
 * @f[ T = \pi \sqrt{\frac{a_t^3}{\mu}} @f]
 *
 * @param h1  Initial orbit altitude above Earth [km]. Must be positive.
 * @param h2  Final   orbit altitude above Earth [km]. Must be positive.
 * @return    HohmannTransfer struct with all computed parameters.
 * @throws std::invalid_argument if either altitude is non-positive.
 */
[[nodiscard]] inline HohmannTransfer compute_hohmann(double h1, double h2)
{
    using namespace constants;

    if (h1 <= 0.0 || h2 <= 0.0)
        throw std::invalid_argument("Orbit altitudes must be positive.");

    const double r1 = R_EARTH + h1;
    const double r2 = R_EARTH + h2;
    const double a  = (r1 + r2) * 0.5;

    const double mu_over_a = MU / a;

    const double vc1 = std::sqrt(MU / r1);
    const double vc2 = std::sqrt(MU / r2);
    const double vt1 = std::sqrt(2.0 * MU / r1 - mu_over_a);
    const double vt2 = std::sqrt(2.0 * MU / r2 - mu_over_a);

    const double dv1 = vt1 - vc1;
    const double dv2 = vc2 - vt2;
    const double tof = PI * a * std::sqrt(a / MU);

    return HohmannTransfer{r1, r2, a, dv1, dv2, tof};
}


/**
 * @brief Circular orbit speed at radius @p r.
 *
 * @f[ v_c = \sqrt{\frac{\mu}{r}} @f]
 *
 * @param r  Orbital radius [km]. Must be positive.
 * @return   Circular velocity [km/s].
 * @throws std::invalid_argument if radius is non-positive.
 */
[[nodiscard]] inline double circular_velocity(double r)
{
    if (r <= 0.0)
        throw std::invalid_argument("Orbital radius must be positive.");
    return std::sqrt(constants::MU / r);
}


/**
 * @brief Keplerian orbital period for semi-major axis @p a.
 *
 * @f[ T = 2\pi \sqrt{\frac{a^3}{\mu}} @f]
 *
 * @param a  Semi-major axis [km]. Must be positive.
 * @return   Orbital period [s].
 * @throws std::invalid_argument if semi-major axis is non-positive.
 */
[[nodiscard]] inline double orbital_period(double a)
{
    if (a <= 0.0)
        throw std::invalid_argument("Semi-major axis must be positive.");
    return 2.0 * constants::PI * a * std::sqrt(a / constants::MU);
}

} // namespace hohmann_het
