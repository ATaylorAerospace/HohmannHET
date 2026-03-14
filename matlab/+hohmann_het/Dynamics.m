% Author: A Taylor | Purpose: Keplerian dynamics and Hohmann transfer | Ref: Vallado/Curtis
classdef Dynamics
% DYNAMICS - Keplerian orbital dynamics and two-impulse Hohmann transfer.
%
% Implements the Hohmann transfer following:
%   Vallado, "Fundamentals of Astrodynamics and Applications," 4th ed.
%   Curtis, "Orbital Mechanics for Engineering Students," 3rd ed.
%
% All lengths in km, velocities in km/s, time in seconds unless noted.
%
% Usage:
%   result = hohmann_het.Dynamics.compute_hohmann(400, 35786)
%   vc     = hohmann_het.Dynamics.circular_velocity(6778.137)
%   T      = hohmann_het.Dynamics.orbital_period(24471.137)
%
% Author: A Taylor | Ref: Vallado/Curtis

    properties (Constant, Access = private)
        %% Physical constants
        MU      = 398600.4418;  % Earth gravitational parameter [km^3/s^2]
        R_EARTH = 6378.137;     % Earth equatorial radius [km]
    end

    methods (Static)

        function result = compute_hohmann(h1, h2)
        % COMPUTE_HOHMANN  Two-impulse Hohmann transfer between circular orbits.
        %
        % Implements Vallado Section 6.3 (impulsive Hohmann transfer).
        %
        % Equations (LaTeX):
        %
        %   r_i = R_E + h_i
        %
        %   a_t = (r_1 + r_2) / 2
        %
        %   v_{c,i} = sqrt(mu / r_i)
        %
        %   v_{t,1} = sqrt(mu * (2/r_1 - 1/a_t)),
        %   v_{t,2} = sqrt(mu * (2/r_2 - 1/a_t))
        %
        %   Delta_v_1 = v_{t,1} - v_{c,1},
        %   Delta_v_2 = v_{c,2} - v_{t,2}
        %
        %   T = pi * sqrt(a_t^3 / mu)
        %
        % Inputs:
        %   h1 - Initial orbit altitude [km], positive finite scalar
        %   h2 - Final   orbit altitude [km], positive finite scalar
        %
        % Output:
        %   result - struct with fields:
        %       r1           - initial orbit radius [km]
        %       r2           - final   orbit radius [km]
        %       a_transfer   - transfer semi-major axis [km]
        %       dv_departure - departure burn delta-V [km/s]
        %       dv_arrival   - arrival  burn delta-V [km/s]
        %       tof          - time of flight [s]
        %       total_dv     - total delta-V [km/s]
        %       tof_hours    - time of flight [h]

            arguments
                h1 (1,1) double {mustBePositive, mustBeFinite}
                h2 (1,1) double {mustBePositive, mustBeFinite}
            end

            MU_     = hohmann_het.Dynamics.MU;
            R_EARTH = hohmann_het.Dynamics.R_EARTH;

            r1 = R_EARTH + h1;
            r2 = R_EARTH + h2;
            a  = (r1 + r2) * 0.5;

            mu_over_a = MU_ / a;

            vc1 = sqrt(MU_ / r1);
            vc2 = sqrt(MU_ / r2);
            vt1 = sqrt(2.0 * MU_ / r1 - mu_over_a);
            vt2 = sqrt(2.0 * MU_ / r2 - mu_over_a);

            dv1 = vt1 - vc1;
            dv2 = vc2 - vt2;
            tof = pi * a * sqrt(a / MU_);

            result.r1           = r1;
            result.r2           = r2;
            result.a_transfer   = a;
            result.dv_departure = dv1;
            result.dv_arrival   = dv2;
            result.tof          = tof;
            result.total_dv     = dv1 + dv2;
            result.tof_hours    = tof / 3600.0;
        end

        function vc = circular_velocity(r)
        % CIRCULAR_VELOCITY  Circular orbit speed at radius r [km].
        %
        %   v_c = sqrt(mu / r)
        %
        % Input:
        %   r  - orbital radius [km], positive finite scalar
        %
        % Output:
        %   vc - circular velocity [km/s]

            arguments
                r (1,1) double {mustBePositive, mustBeFinite}
            end

            vc = sqrt(hohmann_het.Dynamics.MU / r);
        end

        function T = orbital_period(a)
        % ORBITAL_PERIOD  Keplerian orbital period [s].
        %
        %   T = 2*pi * sqrt(a^3 / mu)
        %
        % Input:
        %   a - semi-major axis [km], positive finite scalar
        %
        % Output:
        %   T - orbital period [s]

            arguments
                a (1,1) double {mustBePositive, mustBeFinite}
            end

            T = 2.0 * pi * a * sqrt(a / hohmann_het.Dynamics.MU);
        end

    end % methods (Static)
end % classdef
