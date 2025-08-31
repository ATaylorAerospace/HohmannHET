classdef HohmannTransferOptimized < handle
    properties (Constant, Access = private)
        MU = 398600.4418;        % Earth's gravitational parameter (km^3/s^2)
        R_EARTH = 6378.137;      % Earth's equatorial radius (km)
        INV_3600 = 1/3600;       % Precomputed hours conversion constant
        SPHERE_RESOLUTION = 20;   % Reduced sphere resolution for faster rendering
        ORBIT_POINTS = 100;      % Reduced orbit points for faster plotting
    end
    
    properties (SetAccess = private)
        r1             % Initial orbital radius (km)
        r2             % Final orbital radius (km)
        a_transfer     % Transfer orbit semi-major axis (km)
        delta_v_total  % Total delta-V (km/s)
        transfer_time  % Transfer time (s)
        transfer_hours % Precomputed transfer time in hours
        
        % Precomputed values for display
        initial_altitude
        final_altitude
    end
    
    properties (Access = private)
        delta_v_departure  % Delta-V for departure burn (km/s)
        delta_v_arrival    % Delta-V for arrival burn (km/s)
        
        % Cached trigonometric values for visualization
        cos_theta
        sin_theta
        zeros_array
    end
    
    methods
        function obj = HohmannTransferOptimized(initial_altitude, final_altitude)
            % Simplified validation - faster than validateattributes
            if ~isnumeric(initial_altitude) || ~isscalar(initial_altitude) || ...
               initial_altitude <= 0 || ~isfinite(initial_altitude)
                error('HohmannTransfer:InvalidInput', 'Initial altitude must be a positive finite scalar');
            end
            if ~isnumeric(final_altitude) || ~isscalar(final_altitude) || ...
               final_altitude <= 0 || ~isfinite(final_altitude)
                error('HohmannTransfer:InvalidInput', 'Final altitude must be a positive finite scalar');
            end
            
            % Store for display (avoid repeated calculations)
            obj.initial_altitude = initial_altitude;
            obj.final_altitude = final_altitude;
            
            % Precompute all values at construction
            obj.r1 = obj.R_EARTH + initial_altitude;
            obj.r2 = obj.R_EARTH + final_altitude;
            obj.a_transfer = (obj.r1 + obj.r2) * 0.5;  % Use multiplication instead of division
            
            % Calculate all values in one pass
            obj.computeAllValues();
            
            % Precompute trigonometric values for visualization
            obj.precomputeTrigValues();
        end
        
        function computeAllValues(obj)
            % Vectorized calculations where possible
            mu_over_r = obj.MU ./ [obj.r1, obj.r2];
            v_circular = sqrt(mu_over_r);  % [v1, v2]
            
            % Transfer orbit velocities using optimized calculations
            mu_2_over_r = 2 * mu_over_r;  % [2*MU/r1, 2*MU/r2]
            mu_over_a = obj.MU / obj.a_transfer;
            
            v_transfer = sqrt(mu_2_over_r - mu_over_a);  % [v_transfer_1, v_transfer_2]
            
            % Delta-V calculations
            obj.delta_v_departure = v_transfer(1) - v_circular(1);
            obj.delta_v_arrival = v_circular(2) - v_transfer(2);
            obj.delta_v_total = obj.delta_v_departure + obj.delta_v_arrival;
            
            % Precompute transfer time and hours conversion
            obj.transfer_time = pi * sqrt(obj.a_transfer^3 / obj.MU);
            obj.transfer_hours = obj.transfer_time * obj.INV_3600;
        end
        
        function precomputeTrigValues(obj)
            % Precompute trigonometric values for visualization
            theta = linspace(0, 2*pi, obj.ORBIT_POINTS)';
            obj.cos_theta = cos(theta);
            obj.sin_theta = sin(theta);
            obj.zeros_array = zeros(size(theta));
        end
        
        function printTransferDetails(obj)
            % Single fprintf call with precomputed values
            fprintf(['\nHohmann Transfer Orbit Details:\n', ...
                    'Initial Orbit Altitude: %.2f km\n', ...
                    'Final Orbit Altitude: %.2f km\n', ...
                    'Transfer Orbit Semi-Major Axis: %.2f km\n', ...
                    'Delta-V (Departure Burn): %.2f km/s\n', ...
                    'Delta-V (Arrival Burn): %.2f km/s\n', ...
                    'Total Delta-V: %.2f km/s\n', ...
                    'Transfer Time: %.2f hours\n'], ...
                    obj.initial_altitude, ...
                    obj.final_altitude, ...
                    obj.a_transfer, ...
                    obj.delta_v_departure, ...
                    obj.delta_v_arrival, ...
                    obj.delta_v_total, ...
                    obj.transfer_hours);
        end
        
        function visualizeTransfer(obj)
            % Create figure with optimized settings
            fig = figure('Position', [100, 100, 1200, 800], ...
                        'Renderer', 'opengl');  % Use hardware acceleration
            ax = axes('Parent', fig);
            hold(ax, 'on');
            axis(ax, 'equal');
            
            % Optimized Earth representation with reduced resolution
            [x, y, z] = sphere(obj.SPHERE_RESOLUTION);
            surf(ax, obj.R_EARTH * x, obj.R_EARTH * y, obj.R_EARTH * z, ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'blue', ...
                'DisplayName', 'Earth');
            
            % Plot orbits using precomputed trigonometric values
            plot3(ax, obj.r1 * obj.cos_theta, obj.r1 * obj.sin_theta, obj.zeros_array, ...
                  'g-', 'LineWidth', 2, 'DisplayName', 'Initial Orbit');
            plot3(ax, obj.r2 * obj.cos_theta, obj.r2 * obj.sin_theta, obj.zeros_array, ...
                  'r-', 'LineWidth', 2, 'DisplayName', 'Final Orbit');
            
            % Optimized transfer orbit (ellipse)
            obj.plotTransferOrbit(ax);
            
            % Set visualization properties efficiently
            set(ax, 'XLabel', text(0, 0, 'X (km)'), ...
                   'YLabel', text(0, 0, 'Y (km)'), ...
                   'ZLabel', text(0, 0, 'Z (km)'));
            title(ax, 'Hohmann Transfer Orbit');
            legend(ax, 'show', 'Location', 'best');
            view(ax, [45, 20]);
            grid(ax, 'on');
        end
        
        function plotTransferOrbit(obj, ax)
            % Efficient transfer orbit plotting
            theta_transfer = linspace(0, pi, obj.ORBIT_POINTS/2)';
            
            % Semi-major and semi-minor axes
            a = obj.a_transfer;
            c = abs(obj.r2 - obj.r1) * 0.5;  % Distance from center to focus
            b = sqrt(a^2 - c^2);  % Semi-minor axis
            
            % Ellipse coordinates (centered, then translated)
            x_ellipse = a * cos(theta_transfer);
            y_ellipse = b * sin(theta_transfer);
            
            % Translate to proper position
            x_center = (obj.r1 + obj.r2) * 0.5;
            x_transfer = x_ellipse + x_center - a;
            
            plot3(ax, x_transfer, y_ellipse, zeros(size(theta_transfer)), ...
                  'm--', 'LineWidth', 1.5, 'DisplayName', 'Transfer Orbit');
        end
        
        % Getter methods (inline for speed)
        function val = getDeltaVDeparture(obj)
            val = obj.delta_v_departure;
        end
        
        function val = getDeltaVArrival(obj)
            val = obj.delta_v_arrival;
        end
        
        function val = getTotalDeltaV(obj)
            val = obj.delta_v_total;
        end
        
        function val = getTransferTime(obj)
            val = obj.transfer_time;
        end
        
        function val = getTransferTimeHours(obj)
            val = obj.transfer_hours;
        end
    end
    
    methods (Static)
        function main()
            try
                % More efficient input handling
                initial_alt = HohmannTransferOptimized.getValidInput('Enter initial orbit altitude (km): ');
                final_alt = HohmannTransferOptimized.getValidInput('Enter final orbit altitude (km): ');
                
                transfer = HohmannTransferOptimized(initial_alt, final_alt);
                transfer.printTransferDetails();
                transfer.visualizeTransfer();
            catch ME
                fprintf('Error: %s\n', ME.message);
            end
        end
        
        function value = getValidInput(prompt)
            % Optimized input validation
            while true
                value = input(prompt);
                if isnumeric(value) && isscalar(value) && value > 0 && isfinite(value)
                    break;
                end
                fprintf('Please enter a positive finite number.\n');
            end
        end
    end
end
