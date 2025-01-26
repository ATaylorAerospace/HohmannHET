classdef HohmannTransfer < handle
    properties (Constant, Access = private)
        MU = 398600.4418     % Earth's gravitational parameter (km^3/s^2)
        R_EARTH = 6378.137   % Earth's equatorial radius (km)
    end
    
    properties (SetAccess = private)
        r1             % Initial orbital radius (km)
        r2             % Final orbital radius (km)
        a_transfer     % Transfer orbit semi-major axis (km)
        delta_v_total  % Total delta-V (km/s)
        transfer_time  % Transfer time (s)
    end
    
    properties (Access = private)
        delta_v_departure  % Delta-V for departure burn (km/s)
        delta_v_arrival    % Delta-V for arrival burn (km/s)
    end
    
    methods
        function obj = HohmannTransfer(initial_altitude, final_altitude)
            validateattributes(initial_altitude, {'numeric'}, {'positive', 'scalar', 'finite'}, 'HohmannTransfer', 'initial_altitude');
            validateattributes(final_altitude, {'numeric'}, {'positive', 'scalar', 'finite'}, 'HohmannTransfer', 'final_altitude');
            
            % Precompute all values at construction
            [obj.r1, obj.r2] = deal(obj.R_EARTH + initial_altitude, obj.R_EARTH + final_altitude);
            obj.a_transfer = (obj.r1 + obj.r2) / 2;
            
            % Calculate velocities and delta-Vs
            [obj.delta_v_departure, obj.delta_v_arrival] = obj.computeDeltaVs();
            obj.delta_v_total = obj.delta_v_departure + obj.delta_v_arrival;
            
            % Precompute transfer time
            obj.transfer_time = pi * sqrt(obj.a_transfer^3 / obj.MU);
        end
        
        function [delta_v_dep, delta_v_arr] = computeDeltaVs(obj)
            v1 = sqrt(obj.MU / obj.r1);
            v2 = sqrt(obj.MU / obj.r2);
            
            % Calculate transfer orbit velocities
            v_transfer_1 = sqrt(obj.MU * (2/obj.r1 - 1/obj.a_transfer));
            v_transfer_2 = sqrt(obj.MU * (2/obj.r2 - 1/obj.a_transfer));
            
            delta_v_dep = v_transfer_1 - v1;
            delta_v_arr = v2 - v_transfer_2;
        end
        
        function printTransferDetails(obj)
            fprintf(['\nHohmann Transfer Orbit Details:\n', ...
                    'Initial Orbit Altitude: %.2f km\n', ...
                    'Final Orbit Altitude: %.2f km\n', ...
                    'Transfer Orbit Semi-Major Axis: %.2f km\n', ...
                    'Delta-V (Departure Burn): %.2f km/s\n', ...
                    'Delta-V (Arrival Burn): %.2f km/s\n', ...
                    'Total Delta-V: %.2f km/s\n', ...
                    'Transfer Time: %.2f hours\n'], ...
                    obj.r1 - obj.R_EARTH, ...
                    obj.r2 - obj.R_EARTH, ...
                    obj.a_transfer, ...
                    obj.delta_v_departure, ...
                    obj.delta_v_arrival, ...
                    obj.delta_v_total, ...
                    obj.transfer_time / 3600);
        end
        
        function visualizeTransfer(obj)
            fig = figure('Position', [100, 100, 1200, 800]);
            ax = axes('Parent', fig);
            hold(ax, 'on');
            axis(ax, 'equal');
            
            % Earth representation (more efficient sphere generation)
            [x, y, z] = sphere(30);
            earth = surf(ax, obj.R_EARTH * x, obj.R_EARTH * y, obj.R_EARTH * z, ...
                        'FaceAlpha', 0.1, 'EdgeColor', 'blue', 'DisplayName', 'Earth');
            
            % Generate orbit points efficiently
            theta = linspace(0, 2*pi, 200)';
            cos_theta = cos(theta);
            sin_theta = sin(theta);
            
            % Plot orbits using vectorized operations
            plot3(ax, obj.r1 * cos_theta, obj.r1 * sin_theta, zeros(size(theta)), ...
                  'g-', 'LineWidth', 2, 'DisplayName', 'Initial Orbit');
            plot3(ax, obj.r2 * cos_theta, obj.r2 * sin_theta, zeros(size(theta)), ...
                  'r-', 'LineWidth', 2, 'DisplayName', 'Final Orbit');
            
            % Set visualization properties
            xlabel(ax, 'X (km)');
            ylabel(ax, 'Y (km)');
            zlabel(ax, 'Z (km)');
            title(ax, 'Hohmann Transfer Orbit');
            legend(ax, 'show');
            view(ax, [45, 20]);
            grid(ax, 'on');
        end
    end
    
    methods (Static)
        function main()
            try
                initial_alt = input('Enter initial orbit altitude (km): ');
                final_alt = input('Enter final orbit altitude (km): ');
                
                transfer = HohmannTransfer(initial_alt, final_alt);
                transfer.printTransferDetails();
                transfer.visualizeTransfer();
            catch ME
                fprintf('Error: %s\n', ME.message);
            end
        end
    end
end
