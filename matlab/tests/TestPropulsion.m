% Author: A Taylor | Purpose: matlab.unittest suite for Propulsion | Ref: Vallado/Curtis
classdef TestPropulsion < matlab.unittest.TestCase
% TESTPROPULSION  matlab.unittest test suite for hohmann_het.Propulsion.
%
% Validates HET models against SPT-100 class operating points and
% Tsiolkovsky propellant calculations.
%
% Run with:
%   results = runtests('TestPropulsion');
%   table(results)
%
% Author: A Taylor | Ref: Vallado/Curtis

    properties (Constant)
        TOL      = 1e-6;
        G0       = 9.80665;
        E_CHARGE = 1.602176634e-19;
        M_XE     = 2.180174e-25;

        % SPT-100 reference operating point
        V_D = 300.0;    % V
        P_D = 1350.0;   % W
        ETA = 0.50;     % anode efficiency
    end

    % -----------------------------------------------------------------------
    % Helpers
    % -----------------------------------------------------------------------
    methods (Access = private)

        function ve_ref = ref_ve(tc)
            ve_ref = sqrt(2.0 * tc.ETA * tc.E_CHARGE * tc.V_D / tc.M_XE);
        end

        function isp_ref = ref_isp(tc)
            isp_ref = tc.ref_ve() / tc.G0;
        end

    end

    % -----------------------------------------------------------------------
    % compute_het_state
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'propulsion','het'})

        function testIspParity(tc)
        % Isp matches beam-voltage formula to 1e-6 s.
            s = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            tc.verifyEqual(s.isp, tc.ref_isp(), 'AbsTol', tc.TOL);
        end

        function testExhaustVelocityParity(tc)
        % Exhaust velocity matches sqrt formula to 1e-6 m/s.
            s = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            tc.verifyEqual(s.exhaust_velocity, tc.ref_ve(), 'AbsTol', tc.TOL);
        end

        function testIspReasonableRange(tc)
        % Isp for SPT-100 class should be 1000-2500 s.
            s = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            tc.verifyGreaterThan(s.isp, 1000.0);
            tc.verifyLessThan(s.isp, 2500.0);
        end

        function testThrustPowerRelation(tc)
        % Thrust = 2 * eta * P / ve.
            s  = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            expected_F = 2.0 * tc.ETA * tc.P_D / tc.ref_ve();
            tc.verifyEqual(s.thrust, expected_F, 'AbsTol', tc.TOL);
        end

        function testThrustMassFlowRatio(tc)
        % ve = F / mdot.
            s = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            ve_check = s.thrust / s.mass_flow;
            tc.verifyEqual(ve_check, s.exhaust_velocity, 'AbsTol', tc.TOL);
        end

        function testHigherVoltageHigherIsp(tc)
            s1 = hohmann_het.Propulsion.compute_het_state(300.0, tc.P_D, tc.ETA);
            s2 = hohmann_het.Propulsion.compute_het_state(600.0, tc.P_D, tc.ETA);
            tc.verifyGreaterThan(s2.isp, s1.isp);
        end

        function testStoredEfficiency(tc)
            s = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            tc.verifyEqual(s.anode_efficiency, tc.ETA);
        end

        function testStoredDischargePower(tc)
            s = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            tc.verifyEqual(s.discharge_power, tc.P_D, 'AbsTol', tc.TOL);
        end

    end

    % -----------------------------------------------------------------------
    % Validation errors
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'validation'})

        function testNegativeVoltageThrows(tc)
            tc.verifyError( ...
                @() hohmann_het.Propulsion.compute_het_state(-300, tc.P_D, tc.ETA), ...
                'MATLAB:validators:mustBePositive');
        end

        function testNegativePowerThrows(tc)
            tc.verifyError( ...
                @() hohmann_het.Propulsion.compute_het_state(tc.V_D, -1350, tc.ETA), ...
                'MATLAB:validators:mustBePositive');
        end

        function testEfficiencyZeroThrows(tc)
            tc.verifyError( ...
                @() hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, 0.0), ...
                'MATLAB:validators:mustBeInRange');
        end

        function testEfficiencyOneThrows(tc)
            tc.verifyError( ...
                @() hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, 1.0), ...
                'MATLAB:validators:mustBeInRange');
        end

    end

    % -----------------------------------------------------------------------
    % propellant_mass
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'propulsion','tsiolkovsky'})

        function testPropellantMassParity(tc)
        % Propellant mass matches Tsiolkovsky formula to 1e-6 kg.
            m0  = 1000.0;
            dv  = 3.8526;  % km/s
            isp = tc.ref_isp();
            mp  = hohmann_het.Propulsion.propellant_mass(m0, dv, isp);

            ve       = isp * tc.G0;
            dv_ms    = dv * 1000.0;
            expected = m0 * (1.0 - exp(-dv_ms / ve));
            tc.verifyEqual(mp, expected, 'AbsTol', tc.TOL);
        end

        function testZeroDvZeroProp(tc)
            mp = hohmann_het.Propulsion.propellant_mass(1000.0, 0.0, tc.ref_isp());
            tc.verifyEqual(mp, 0.0, 'AbsTol', tc.TOL);
        end

        function testHigherIspLowerPropMass(tc)
            mp1 = hohmann_het.Propulsion.propellant_mass(1000.0, 3.0, 1500.0);
            mp2 = hohmann_het.Propulsion.propellant_mass(1000.0, 3.0, 3000.0);
            tc.verifyLessThan(mp2, mp1);
        end

        function testPropMassLessThanInitial(tc)
            mp = hohmann_het.Propulsion.propellant_mass(1000.0, 3.8526, tc.ref_isp());
            tc.verifyLessThan(mp, 1000.0);
        end

    end

    % -----------------------------------------------------------------------
    % burn_time
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'propulsion'})

        function testBurnTimePositive(tc)
            s  = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            tb = hohmann_het.Propulsion.burn_time(s.thrust, s.mass_flow, 3.8526, 1000.0);
            tc.verifyGreaterThan(tb, 0.0);
        end

        function testLargerDvLongerBurn(tc)
            s   = hohmann_het.Propulsion.compute_het_state(tc.V_D, tc.P_D, tc.ETA);
            tb1 = hohmann_het.Propulsion.burn_time(s.thrust, s.mass_flow, 1.0, 1000.0);
            tb2 = hohmann_het.Propulsion.burn_time(s.thrust, s.mass_flow, 3.0, 1000.0);
            tc.verifyGreaterThan(tb2, tb1);
        end

    end

end
