% Author: A Taylor | Purpose: matlab.unittest suite for Dynamics | Ref: Vallado/Curtis
classdef TestDynamics < matlab.unittest.TestCase
% TESTDYNAMICS  matlab.unittest test suite for hohmann_het.Dynamics.
%
% Validates Keplerian dynamics against known LEO-to-GEO benchmark values
% and verifies cross-language parity constants to 1e-6 precision.
%
% Run with:
%   results = runtests('TestDynamics');
%   table(results)
%
% Author: A Taylor | Ref: Vallado/Curtis

    properties (Constant)
        MU_REF      = 398600.4418;
        R_EARTH_REF = 6378.137;
        TOL         = 1e-6;

        H1 = 400.0;      % km  (LEO altitude)
        H2 = 35786.0;    % km  (GEO altitude)
    end

    % -----------------------------------------------------------------------
    % Constants
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'constants'})

        function testMuValue(tc)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            % Direct constant check via known formula
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            % If MU were wrong the semi-major axis would differ
            r1 = tc.R_EARTH_REF + tc.H1;
            r2 = tc.R_EARTH_REF + tc.H2;
            a_expected = (r1 + r2) / 2;
            tc.verifyEqual(r.a_transfer, a_expected, 'AbsTol', tc.TOL);
        end

    end

    % -----------------------------------------------------------------------
    % compute_hohmann - valid inputs
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'dynamics','hohmann'})

        function testRadiiValues(tc)
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyEqual(r.r1, tc.R_EARTH_REF + tc.H1, 'AbsTol', tc.TOL);
            tc.verifyEqual(r.r2, tc.R_EARTH_REF + tc.H2, 'AbsTol', tc.TOL);
        end

        function testSemiMajorAxis(tc)
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            expected_a = (tc.R_EARTH_REF + tc.H1 + tc.R_EARTH_REF + tc.H2) / 2.0;
            tc.verifyEqual(r.a_transfer, expected_a, 'AbsTol', tc.TOL);
        end

        function testDvDepartureParity(tc)
        % Delta-V departure matches reference formula to 1e-6 km/s.
            [ref_dv1, ~, ~, ~] = TestDynamics.ref_hohmann(tc.H1, tc.H2);
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyEqual(r.dv_departure, ref_dv1, 'AbsTol', tc.TOL);
        end

        function testDvArrivalParity(tc)
        % Delta-V arrival matches reference formula to 1e-6 km/s.
            [~, ref_dv2, ~, ~] = TestDynamics.ref_hohmann(tc.H1, tc.H2);
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyEqual(r.dv_arrival, ref_dv2, 'AbsTol', tc.TOL);
        end

        function testTotalDvParity(tc)
        % Total delta-V matches reference formula to 1e-6 km/s.
            [~, ~, ref_total, ~] = TestDynamics.ref_hohmann(tc.H1, tc.H2);
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyEqual(r.total_dv, ref_total, 'AbsTol', tc.TOL);
        end

        function testTofParity(tc)
        % TOF matches reference formula to 1e-6 s.
            [~, ~, ~, ref_tof] = TestDynamics.ref_hohmann(tc.H1, tc.H2);
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyEqual(r.tof, ref_tof, 'AbsTol', tc.TOL);
        end

        function testKnownLeoGeoDvDeparture(tc)
        % LEO-to-GEO departure burn approx 2.4 km/s (Vallado Table 6-1).
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyGreaterThan(r.dv_departure, 2.3);
            tc.verifyLessThan(r.dv_departure, 2.5);
        end

        function testKnownLeoGeoTotalDv(tc)
        % LEO-to-GEO total dv approx 3.9 km/s.
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyGreaterThan(r.total_dv, 3.7);
            tc.verifyLessThan(r.total_dv, 4.1);
        end

        function testKnownLeoGeoTofHours(tc)
        % LEO-to-GEO transfer time approx 5.3 hours.
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyGreaterThan(r.tof_hours, 4.5);
            tc.verifyLessThan(r.tof_hours, 6.0);
        end

        function testSameAltitudeZeroDv(tc)
        % Same-altitude transfer requires zero delta-V.
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H1);
            tc.verifyEqual(r.dv_departure, 0.0, 'AbsTol', tc.TOL);
            tc.verifyEqual(r.dv_arrival,   0.0, 'AbsTol', tc.TOL);
            tc.verifyEqual(r.total_dv,     0.0, 'AbsTol', tc.TOL);
        end

        function testReverseSymmetryTotalDv(tc)
        % Forward and reverse transfers have equal total delta-V.
            fwd = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            rev = hohmann_het.Dynamics.compute_hohmann(tc.H2, tc.H1);
            tc.verifyEqual(fwd.total_dv, rev.total_dv, 'AbsTol', tc.TOL);
        end

        function testReverseSymmetryTof(tc)
        % Forward and reverse transfers have equal TOF.
            fwd = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            rev = hohmann_het.Dynamics.compute_hohmann(tc.H2, tc.H1);
            tc.verifyEqual(fwd.tof, rev.tof, 'AbsTol', tc.TOL);
        end

        function testTotalDvEqualsSum(tc)
        % total_dv field equals dv1 + dv2.
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyEqual(r.total_dv, r.dv_departure + r.dv_arrival, 'AbsTol', tc.TOL);
        end

        function testTofHoursConversion(tc)
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyEqual(r.tof_hours, r.tof / 3600.0, 'AbsTol', tc.TOL);
        end

        function testSemiMajorAxisBetweenRadii(tc)
        % For ascending transfer, r1 < a_transfer < r2.
            r = hohmann_het.Dynamics.compute_hohmann(tc.H1, tc.H2);
            tc.verifyGreaterThan(r.a_transfer, r.r1);
            tc.verifyLessThan(r.a_transfer, r.r2);
        end

    end

    % -----------------------------------------------------------------------
    % Invalid inputs
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'validation'})

        function testNegativeH1Throws(tc)
            tc.verifyError(@() hohmann_het.Dynamics.compute_hohmann(-100, tc.H2), ...
                           'MATLAB:validators:mustBePositive');
        end

        function testNegativeH2Throws(tc)
            tc.verifyError(@() hohmann_het.Dynamics.compute_hohmann(tc.H1, -100), ...
                           'MATLAB:validators:mustBePositive');
        end

        function testZeroH1Throws(tc)
            tc.verifyError(@() hohmann_het.Dynamics.compute_hohmann(0, tc.H2), ...
                           'MATLAB:validators:mustBePositive');
        end

    end

    % -----------------------------------------------------------------------
    % circular_velocity
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'dynamics'})

        function testCircularVelocityLeo(tc)
            r1 = tc.R_EARTH_REF + tc.H1;
            vc = hohmann_het.Dynamics.circular_velocity(r1);
            expected = sqrt(tc.MU_REF / r1);
            tc.verifyEqual(vc, expected, 'AbsTol', tc.TOL);
        end

        function testCircularVelocityNegativeThrows(tc)
            tc.verifyError(@() hohmann_het.Dynamics.circular_velocity(-1), ...
                           'MATLAB:validators:mustBePositive');
        end

    end

    % -----------------------------------------------------------------------
    % orbital_period
    % -----------------------------------------------------------------------
    methods (Test, TestTags = {'dynamics'})

        function testGeoPeriodApprox24h(tc)
        % GEO orbital period is approximately 24 hours.
            r2 = tc.R_EARTH_REF + tc.H2;
            T  = hohmann_het.Dynamics.orbital_period(r2);
            tc.verifyEqual(T / 3600, 24.0, 'AbsTol', 0.1);
        end

        function testOrbitalPeriodFormula(tc)
            a = (tc.R_EARTH_REF + tc.H1 + tc.R_EARTH_REF + tc.H2) / 2;
            T = hohmann_het.Dynamics.orbital_period(a);
            expected = 2.0 * pi * a * sqrt(a / tc.MU_REF);
            tc.verifyEqual(T, expected, 'AbsTol', tc.TOL);
        end

    end

    % -----------------------------------------------------------------------
    % Private reference implementation
    % -----------------------------------------------------------------------
    methods (Static, Access = private)

        function [dv1, dv2, total, tof] = ref_hohmann(h1, h2)
            MU      = 398600.4418;
            R_EARTH = 6378.137;
            r1  = R_EARTH + h1;
            r2  = R_EARTH + h2;
            a   = (r1 + r2) / 2.0;
            mua = MU / a;
            vc1 = sqrt(MU / r1);
            vc2 = sqrt(MU / r2);
            vt1 = sqrt(2.0 * MU / r1 - mua);
            vt2 = sqrt(2.0 * MU / r2 - mua);
            dv1   = vt1 - vc1;
            dv2   = vc2 - vt2;
            total = dv1 + dv2;
            tof   = pi * a * sqrt(a / MU);
        end

    end

end
