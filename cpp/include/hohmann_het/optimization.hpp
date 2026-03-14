// Author: A Taylor | Purpose: Mission optimization for HET transfers | Ref: Vallado/Curtis
#pragma once

/**
 * @file optimization.hpp
 * @brief Mission optimization solvers for HET orbital transfers.
 *
 * Implements golden-section search for minimizing a composite
 * propellant-mass / burn-time objective over specific impulse.
 *
 * The golden-section algorithm is numerically identical to the Python
 * and MATLAB implementations, guaranteeing cross-language result parity
 * within floating-point tolerance (1e-6).
 *
 * @author A Taylor
 * @ref Vallado/Curtis
 */

#include <cmath>
#include <functional>
#include <stdexcept>

#include "dynamics.hpp"
#include "propulsion.hpp"

namespace hohmann_het {

// ---------------------------------------------------------------------------
// Result struct
// ---------------------------------------------------------------------------
/**
 * @brief Result of a specific-impulse optimization.
 */
struct OptimizationResult {
    double optimal_isp;      ///< Optimal specific impulse [s]
    double propellant_mass;  ///< Propellant mass at the optimum [kg]
    double burn_time;        ///< Burn time at the optimum [s]
    double objective_value;  ///< Composite objective value at the optimum
};


// ---------------------------------------------------------------------------
// Core algorithm
// ---------------------------------------------------------------------------
/**
 * @brief Golden-section search for the minimum of a unimodal scalar function.
 *
 * Converges to within @p tol of the true minimum location.
 *
 * The interval is reduced at each step by the golden ratio:
 * @f[ \varphi = \frac{1 + \sqrt{5}}{2} \approx 1.618 @f]
 *
 * @param func  Scalar function to minimize.
 * @param a     Left  bound of search interval.
 * @param b     Right bound of search interval.
 * @param tol   Convergence tolerance on interval width.
 * @return      Approximate location of the minimum.
 */
inline double golden_section_minimize(
    std::function<double(double)> func,
    double a,
    double b,
    double tol = 1e-9)
{
    const double gr = (std::sqrt(5.0) + 1.0) / 2.0;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;

    while (std::abs(b - a) > tol) {
        if (func(c) < func(d))
            b = d;
        else
            a = c;
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    }
    return (a + b) / 2.0;
}


// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------
/**
 * @brief Optimize specific impulse to minimize a composite propellant-time objective.
 *
 * @par Objective function
 * @f[
 *   J(I_{sp}) = \frac{m_p(I_{sp})}{m_0}
 *             + \lambda \frac{T_b(I_{sp})}{T_{\rm ref}}
 * @f]
 *
 * Propellant mass (Tsiolkovsky):
 * @f[ m_p = m_0 \left(1 - e^{-\Delta v / (g_0 I_{sp})}\right) @f]
 *
 * Power-limited burn time:
 * @f[ T_b = \frac{m_0 \,\Delta v \, v_e}{2 \eta_a P_d}, \quad v_e = g_0 I_{sp} @f]
 *
 * Minimized via golden-section search over
 * @f$I_{sp} \in [I_{sp,\min},\ I_{sp,\max}]@f$.
 *
 * @param m_initial         Spacecraft wet mass [kg]. Must be positive.
 * @param delta_v           Required delta-V [m/s]. Non-negative.
 * @param discharge_power   Available thruster power [W]. Must be positive.
 * @param anode_efficiency  HET anode efficiency @f$\eta_a \in (0,1)@f$.
 * @param isp_min           Lower bound of Isp search [s]. Default 500.
 * @param isp_max           Upper bound of Isp search [s]. Default 5000.
 * @param lambda_weight     Burn-time penalty weight. Default 1e-4.
 * @return                  OptimizationResult with optimal Isp and mission parameters.
 * @throws std::invalid_argument on invalid inputs.
 */
[[nodiscard]] inline OptimizationResult optimize_isp(
    double m_initial,
    double delta_v,
    double discharge_power,
    double anode_efficiency,
    double isp_min        = 500.0,
    double isp_max        = 5000.0,
    double lambda_weight  = 1e-4)
{
    if (m_initial <= 0.0 || delta_v < 0.0 || discharge_power <= 0.0)
        throw std::invalid_argument("Mass, delta-V, and power must be non-negative/positive.");
    if (anode_efficiency <= 0.0 || anode_efficiency >= 1.0)
        throw std::invalid_argument("Anode efficiency must be in (0, 1).");
    if (isp_min <= 0.0 || isp_max <= isp_min)
        throw std::invalid_argument("Isp bounds must satisfy 0 < isp_min < isp_max.");

    const double g0    = constants::G0;
    const double ve_ref = isp_min * g0;
    const double T_ref  = m_initial * delta_v * ve_ref
                          / (2.0 * anode_efficiency * discharge_power);

    auto objective = [&](double isp_s) -> double {
        const double ve      = isp_s * g0;
        const double mp_frac = 1.0 - std::exp(-delta_v / ve);
        const double t_burn  = m_initial * delta_v * ve
                               / (2.0 * anode_efficiency * discharge_power);
        return mp_frac + lambda_weight * (t_burn / T_ref);
    };

    const double isp_opt = golden_section_minimize(objective, isp_min, isp_max);
    const double ve_opt  = isp_opt * g0;
    const double mp_opt  = m_initial * (1.0 - std::exp(-delta_v / ve_opt));
    const double tb_opt  = m_initial * delta_v * ve_opt
                           / (2.0 * anode_efficiency * discharge_power);
    const double J_opt   = objective(isp_opt);

    return OptimizationResult{isp_opt, mp_opt, tb_opt, J_opt};
}

} // namespace hohmann_het
