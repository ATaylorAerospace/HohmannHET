% Author: A Taylor | Purpose: matlab.unittest suite for Optimization | Ref: Vallado/Curtis
classdef TestOptimization < matlab.unittest.TestCase
% TESTOPTIMIZATION  matlab.unittest test suite for hohmann_het.Optimization.
%
% Validates golden-section search and mission optimization against
% expected analytical behavior.
%
% Run with:
%   results = runtests('TestOptimization');
%   table(results)
%
% Author: A Taylor | Ref: Vallado/Curtis

    properties (Constant)
        TOL = 1e-6;
        M0  = 1000.0;    % kg
        DV  = 3852.6;    % m/s  (LEO-GEO delta-V)
        P   = 5000.0;    % W
        ETA = 0.55;
    end

    % -----------------------------------------------------------------------
    % golden_section_minimize
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'algorithm'})

        function testQuadraticMinimum(tc)
        % Find minimum of f(x) = (x-3)^2 on [0, 10].
            f = @(x) (x - 3.0)^2;
            x_min = hohmann_het.Optimization.golden_section_minimize(f, 0.0, 10.0);
            tc.verifyEqual(x_min, 3.0, 'AbsTol', 1e-6);
        end

        function testOffsetQuadratic(tc)
            f = @(x) (x - 1.5)^2;
            x_min = hohmann_het.Optimization.golden_section_minimize(f, -5.0, 5.0);
            tc.verifyEqual(x_min, 1.5, 'AbsTol', 1e-6);
        end

        function testTighterTolerance(tc)
            f = @(x) (x - 2.718281828)^2;
            x_min = hohmann_het.Optimization.golden_section_minimize(f, 0.0, 5.0, 1e-12);
            tc.verifyEqual(x_min, 2.718281828, 'AbsTol', 1e-9);
        end

    end

    % -----------------------------------------------------------------------
    % optimize_isp
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'optimization'})

        function testIspWithinBounds(tc)
            r = hohmann_het.Optimization.optimize_isp( ...
                tc.M0, tc.DV, tc.P, tc.ETA, 500.0, 5000.0);
            tc.verifyGreaterThanOrEqual(r.optimal_isp, 500.0);
            tc.verifyLessThanOrEqual(r.optimal_isp, 5000.0);
        end

        function testPropellantMassPositive(tc)
            r = hohmann_het.Optimization.optimize_isp(tc.M0, tc.DV, tc.P, tc.ETA);
            tc.verifyGreaterThan(r.propellant_mass, 0.0);
        end

        function testPropellantMassLessThanInitial(tc)
            r = hohmann_het.Optimization.optimize_isp(tc.M0, tc.DV, tc.P, tc.ETA);
            tc.verifyLessThan(r.propellant_mass, tc.M0);
        end

        function testBurnTimePositive(tc)
            r = hohmann_het.Optimization.optimize_isp(tc.M0, tc.DV, tc.P, tc.ETA);
            tc.verifyGreaterThan(r.burn_time, 0.0);
        end

        function testObjectiveValueFinite(tc)
            r = hohmann_het.Optimization.optimize_isp(tc.M0, tc.DV, tc.P, tc.ETA);
            tc.verifyTrue(isfinite(r.objective_value));
        end

        function testHigherLambdaLowersIsp(tc)
        % Heavier burn-time penalty pushes optimizer toward lower Isp.
            r_low  = hohmann_het.Optimization.optimize_isp( ...
                tc.M0, tc.DV, tc.P, tc.ETA, 500, 5000, 1e-6);
            r_high = hohmann_het.Optimization.optimize_isp( ...
                tc.M0, tc.DV, tc.P, tc.ETA, 500, 5000, 1e-1);
            tc.verifyLessThanOrEqual(r_high.optimal_isp, r_low.optimal_isp);
        end

        function testBadBoundsThrows(tc)
            tc.verifyError( ...
                @() hohmann_het.Optimization.optimize_isp( ...
                    tc.M0, tc.DV, tc.P, tc.ETA, 3000.0, 1000.0), ...
                'hohmann_het:Optimization:badBounds');
        end

    end

    % -----------------------------------------------------------------------
    % min_propellant_transfer
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'optimization','integration'})

        function testEndToEndLeoGeo(tc)
        % Wrapper produces positive propellant mass for LEO-GEO mission.
            r = hohmann_het.Optimization.min_propellant_transfer( ...
                400.0, 35786.0, tc.M0, tc.P, tc.ETA);
            tc.verifyGreaterThan(r.propellant_mass, 0.0);
            tc.verifyGreaterThan(r.optimal_isp, 0.0);
        end

        function testConsistencyWithOptimizeIsp(tc)
        % Wrapper result matches direct optimize_isp call.
            transfer = hohmann_het.Dynamics.compute_hohmann(400.0, 35786.0);
            dv_ms    = transfer.total_dv * 1000.0;

            r_direct = hohmann_het.Optimization.optimize_isp( ...
                tc.M0, dv_ms, tc.P, tc.ETA);
            r_wrap   = hohmann_het.Optimization.min_propellant_transfer( ...
                400.0, 35786.0, tc.M0, tc.P, tc.ETA);

            tc.verifyEqual(r_direct.optimal_isp,    r_wrap.optimal_isp,    'AbsTol', tc.TOL);
            tc.verifyEqual(r_direct.propellant_mass, r_wrap.propellant_mass, 'AbsTol', tc.TOL);
        end

    end

end
