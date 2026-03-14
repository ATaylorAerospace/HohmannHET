% Author: A Taylor | Purpose: Hall Effect Thruster propulsion models | Ref: Vallado/Curtis
classdef Propulsion
% PROPULSION - High-fidelity Hall Effect Thruster (HET) models.
%
% Models anode efficiency, beam voltage, exhaust velocity, specific impulse,
% mass flow rate, and propellant budget.
%
% Reference:
%   Goebel & Katz, "Fundamentals of Electric Propulsion," JPL/Wiley, 2008.
%   Vallado, "Fundamentals of Astrodynamics and Applications," 4th ed.
%
% Usage:
%   state = hohmann_het.Propulsion.compute_het_state(300, 1350, 0.50)
%   mp    = hohmann_het.Propulsion.propellant_mass(1000, 3.8526, 1514)
%
% Author: A Taylor | Ref: Vallado/Curtis

    properties (Constant, Access = private)
        %% Physical constants
        G0       = 9.80665;          % Standard gravity [m/s^2]
        E_CHARGE = 1.602176634e-19;  % Elementary charge [C]
        M_XE     = 2.180174e-25;     % Xenon atom mass [kg]  (131.293 u)
    end

    methods (Static)

        function state = compute_het_state(discharge_voltage, discharge_power, anode_efficiency)
        % COMPUTE_HET_STATE  Compute HET operating point from discharge parameters.
        %
        % Equations (LaTeX):
        %
        %   v_e = sqrt(2 * eta_a * e * V_d / m_Xe)   [beam exhaust velocity, m/s]
        %
        %   I_sp = v_e / g_0                           [specific impulse, s]
        %
        %   F = 2 * eta_a * P_d / v_e                 [thrust, N]
        %
        %   m_dot = F / v_e                            [mass flow rate, kg/s]
        %
        % Inputs:
        %   discharge_voltage  - anode discharge voltage [V],  positive scalar
        %   discharge_power    - input discharge power   [W],  positive scalar
        %   anode_efficiency   - anode efficiency eta_a in (0,1), scalar
        %
        % Output:
        %   state - struct with fields:
        %       isp              [s]
        %       exhaust_velocity [m/s]
        %       thrust           [N]
        %       mass_flow        [kg/s]
        %       anode_efficiency [-]
        %       discharge_power  [W]

            arguments
                discharge_voltage  (1,1) double {mustBePositive, mustBeFinite}
                discharge_power    (1,1) double {mustBePositive, mustBeFinite}
                anode_efficiency   (1,1) double {mustBeInRange(anode_efficiency, 0, 1, "exclusive", "exclusive")}
            end

            G0       = hohmann_het.Propulsion.G0;
            E_CHARGE = hohmann_het.Propulsion.E_CHARGE;
            M_XE     = hohmann_het.Propulsion.M_XE;

            ve   = sqrt(2.0 * anode_efficiency * E_CHARGE * discharge_voltage / M_XE);
            isp  = ve / G0;
            F    = 2.0 * anode_efficiency * discharge_power / ve;
            mdot = F / ve;

            state.isp              = isp;
            state.exhaust_velocity = ve;
            state.thrust           = F;
            state.mass_flow        = mdot;
            state.anode_efficiency = anode_efficiency;
            state.discharge_power  = discharge_power;
        end

        function mp = propellant_mass(m_initial, delta_v_km_s, isp)
        % PROPELLANT_MASS  Tsiolkovsky rocket equation.
        %
        %   m_p = m_0 * (1 - exp(-Delta_v / (g_0 * I_sp)))
        %
        % Inputs:
        %   m_initial    - initial spacecraft mass [kg],  positive scalar
        %   delta_v_km_s - required delta-V [km/s],       non-negative scalar
        %   isp          - specific impulse [s],           positive scalar
        %
        % Output:
        %   mp - propellant mass [kg]

            arguments
                m_initial    (1,1) double {mustBePositive,    mustBeFinite}
                delta_v_km_s (1,1) double {mustBeNonnegative, mustBeFinite}
                isp          (1,1) double {mustBePositive,    mustBeFinite}
            end

            G0 = hohmann_het.Propulsion.G0;
            ve = isp * G0;                      % m/s
            dv = delta_v_km_s * 1000.0;         % km/s -> m/s
            mp = m_initial * (1.0 - exp(-dv / ve));
        end

        function tb = burn_time(thrust, mass_flow, delta_v_km_s, m_initial)
        % BURN_TIME  Constant-thrust burn time for a given delta-V.
        %
        %   T_b = (m_0 / m_dot) * (1 - exp(-Delta_v / v_e))
        %
        % Inputs:
        %   thrust       - thruster force [N],      positive scalar
        %   mass_flow    - mass flow rate [kg/s],   positive scalar
        %   delta_v_km_s - required delta-V [km/s], non-negative scalar
        %   m_initial    - initial mass [kg],        positive scalar
        %
        % Output:
        %   tb - burn time [s]

            arguments
                thrust       (1,1) double {mustBePositive,    mustBeFinite}
                mass_flow    (1,1) double {mustBePositive,    mustBeFinite}
                delta_v_km_s (1,1) double {mustBeNonnegative, mustBeFinite}
                m_initial    (1,1) double {mustBePositive,    mustBeFinite}
            end

            ve = thrust / mass_flow;                % m/s
            dv = delta_v_km_s * 1000.0;             % km/s -> m/s
            tb = (m_initial / mass_flow) * (1.0 - exp(-dv / ve));
        end

    end % methods (Static)
end % classdef
