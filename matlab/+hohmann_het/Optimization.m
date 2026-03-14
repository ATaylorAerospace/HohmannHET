% Author: A Taylor | Purpose: Mission optimization solvers for HET transfers | Ref: Vallado/Curtis
classdef Optimization
% OPTIMIZATION - Mission optimization solvers for HET orbital transfers.
%
% Implements golden-section search for minimizing a composite
% propellant-mass / burn-time objective over specific impulse.
%
% The golden-section algorithm is numerically identical to the Python
% and C++ implementations, guaranteeing cross-language result parity
% within floating-point tolerance (1e-6).
%
% Usage:
%   r = hohmann_het.Optimization.optimize_isp(1000, 3852.6, 5000, 0.55)
%   r = hohmann_het.Optimization.min_propellant_transfer(400, 35786, 1000, 5000, 0.55)
%
% Author: A Taylor | Ref: Vallado/Curtis

    properties (Constant, Access = private)
        G0 = 9.80665;  % Standard gravity [m/s^2]
    end

    methods (Static)

        function x_min = golden_section_minimize(func, a, b, tol)
        % GOLDEN_SECTION_MINIMIZE  Minimize a unimodal scalar function on [a, b].
        %
        % Converges to within tol of the true minimum location.
        %
        %   phi = (1 + sqrt(5)) / 2   (golden ratio)
        %
        % Inputs:
        %   func - function handle @(x) scalar
        %   a    - left  bound of search interval
        %   b    - right bound of search interval
        %   tol  - convergence tolerance (default 1e-9)
        %
        % Output:
        %   x_min - approximate location of the minimum

            arguments
                func (1,1) {mustBeA(func, "function_handle")}
                a    (1,1) double {mustBeFinite}
                b    (1,1) double {mustBeFinite}
                tol  (1,1) double {mustBePositive} = 1e-9
            end

            gr = (sqrt(5.0) + 1.0) / 2.0;  % golden ratio
            c  = b - (b - a) / gr;
            d  = a + (b - a) / gr;

            while abs(b - a) > tol
                if func(c) < func(d)
                    b = d;
                else
                    a = c;
                end
                c = b - (b - a) / gr;
                d = a + (b - a) / gr;
            end

            x_min = (a + b) / 2.0;
        end

        function result = optimize_isp(m_initial, delta_v_m_s, discharge_power, ...
                                        anode_efficiency, isp_min, isp_max, lambda_weight)
        % OPTIMIZE_ISP  Find specific impulse minimizing composite propellant-time objective.
        %
        % Objective function:
        %
        %   J(I_sp) = m_p(I_sp)/m_0  +  lambda * T_b(I_sp)/T_ref
        %
        % Propellant mass (Tsiolkovsky):
        %
        %   m_p = m_0 * (1 - exp(-Delta_v / (g_0 * I_sp)))
        %
        % Power-limited burn time:
        %
        %   T_b = m_0 * Delta_v * v_e / (2 * eta_a * P_d),  v_e = g_0 * I_sp
        %
        % Minimized via golden-section search over I_sp in [isp_min, isp_max].
        %
        % Inputs:
        %   m_initial        - spacecraft wet mass [kg],   positive scalar
        %   delta_v_m_s      - required delta-V [m/s],     non-negative scalar
        %   discharge_power  - thruster power [W],          positive scalar
        %   anode_efficiency - eta_a in (0,1),              scalar
        %   isp_min          - lower Isp bound [s]          (default 500)
        %   isp_max          - upper Isp bound [s]          (default 5000)
        %   lambda_weight    - burn-time penalty weight     (default 1e-4)
        %
        % Output:
        %   result - struct with fields:
        %       optimal_isp      [s]
        %       propellant_mass  [kg]
        %       burn_time        [s]
        %       objective_value  [-]

            arguments
                m_initial        (1,1) double {mustBePositive,    mustBeFinite}
                delta_v_m_s      (1,1) double {mustBeNonnegative, mustBeFinite}
                discharge_power  (1,1) double {mustBePositive,    mustBeFinite}
                anode_efficiency (1,1) double {mustBeInRange(anode_efficiency, 0, 1, "exclusive", "exclusive")}
                isp_min          (1,1) double {mustBePositive,    mustBeFinite} = 500.0
                isp_max          (1,1) double {mustBePositive,    mustBeFinite} = 5000.0
                lambda_weight    (1,1) double {mustBePositive,    mustBeFinite} = 1e-4
            end

            if isp_max <= isp_min
                error('hohmann_het:Optimization:badBounds', ...
                      'isp_max must be greater than isp_min.');
            end

            G0    = hohmann_het.Optimization.G0;
            m0    = m_initial;
            dv    = delta_v_m_s;
            P     = discharge_power;
            eta_a = anode_efficiency;

            % Reference burn time at minimum Isp (for normalisation)
            ve_ref = isp_min * G0;
            T_ref  = m0 * dv * ve_ref / (2.0 * eta_a * P);

            obj = @(isp_s) ...
                (1.0 - exp(-dv / (isp_s * G0))) + ...
                lambda_weight * (m0 * dv * (isp_s * G0) / (2.0 * eta_a * P)) / T_ref;

            isp_opt = hohmann_het.Optimization.golden_section_minimize( ...
                obj, isp_min, isp_max);

            ve_opt = isp_opt * G0;
            mp_opt = m0 * (1.0 - exp(-dv / ve_opt));
            tb_opt = m0 * dv * ve_opt / (2.0 * eta_a * P);
            J_opt  = obj(isp_opt);

            result.optimal_isp     = isp_opt;
            result.propellant_mass = mp_opt;
            result.burn_time       = tb_opt;
            result.objective_value = J_opt;
        end

        function result = min_propellant_transfer(h1, h2, m_initial, ...
                                                   discharge_power, anode_efficiency, ...
                                                   isp_min, isp_max)
        % MIN_PROPELLANT_TRANSFER  Compute Hohmann dv then optimize Isp.
        %
        % Convenience wrapper combining Dynamics.compute_hohmann and optimize_isp.
        %
        % Inputs:
        %   h1, h2           - initial/final orbit altitudes [km]
        %   m_initial        - spacecraft wet mass [kg]
        %   discharge_power  - thruster power [W]
        %   anode_efficiency - eta_a in (0,1)
        %   isp_min          - lower Isp bound [s]  (default 500)
        %   isp_max          - upper Isp bound [s]  (default 5000)
        %
        % Output:
        %   result - same struct as optimize_isp

            arguments
                h1               (1,1) double {mustBePositive, mustBeFinite}
                h2               (1,1) double {mustBePositive, mustBeFinite}
                m_initial        (1,1) double {mustBePositive, mustBeFinite}
                discharge_power  (1,1) double {mustBePositive, mustBeFinite}
                anode_efficiency (1,1) double {mustBeInRange(anode_efficiency, 0, 1, "exclusive", "exclusive")}
                isp_min          (1,1) double {mustBePositive, mustBeFinite} = 500.0
                isp_max          (1,1) double {mustBePositive, mustBeFinite} = 5000.0
            end

            transfer = hohmann_het.Dynamics.compute_hohmann(h1, h2);
            dv_m_s   = transfer.total_dv * 1000.0;   % km/s -> m/s

            result = hohmann_het.Optimization.optimize_isp( ...
                m_initial, dv_m_s, discharge_power, anode_efficiency, ...
                isp_min, isp_max);
        end

    end % methods (Static)
end % classdef
