classdef HohmannTransfer
    % Calculation Constants
    properties (Constant)
        MU = 398600.4418;     % Earth's gravitational parameter (km^3/s^2)
        R_EARTH = 6378.137;   % Earth's equatorial radius (km)
    end
    
    % Private properties to store transfer orbit parameters
    properties (Access = private)
        r1             % Initial orbital radius (km)
        r2             % Final orbital radius (km)
        a_transfer     % Transfer orbit semi-major axis (km)
        v1             % Initial orbital velocity (km/s)
        delta_v_departure  % Delta-V for departure burn (km/s)
        delta_v_arrival    % Delta-V for arrival burn (km/s)
    end
    
    methods
        % Constructor
        function obj = HohmannTransfer(initial_altitude, final_altitude)
            % Validate inputs
            if initial_altitude <= 0 || final_altitude <= 0
                error('Altitudes must be positive');
            end
            
            % Calculate orbital radii
            obj.r1 = obj.R_EARTH + initial_altitude;
            obj.r2 = obj.R_EARTH + final_altitude;
            
            % Calculate transfer orbit semi-major axis
            obj.a_transfer = (obj.r1 + obj.r2) / 2;
            
            % Calculate initial orbital velocity
            obj.v1 = sqrt(obj.MU / obj.r1);
            
            % Compute delta-V for transfer
            obj.delta_v_departure = obj.calculateDeltaVDeparture();
            obj.delta_v_arrival = obj.calculateDeltaVArrival();
        end
        
        % Calculate delta-V for departure burn
        function delta_v = calculateDeltaVDeparture(obj)
            v_transfer_1 = sqrt(obj.MU * (2/obj.r1 - 1/obj.a_transfer));
            delta_v = v_transfer_1 - obj.v1;
        end
        
        % Calculate delta-V for arrival burn
        function delta_v = calculateDeltaVArrival(obj)
            % Velocity at final orbit
            v2 = sqrt(obj.MU / obj.r2);
            
            % Velocity at transfer orbit apoapsis
            v_transfer_2 = sqrt(obj.MU * (2/obj.r2 - 1/obj.a_transfer));
            
            delta_v = v2 - v_transfer_2;
        end
        
        % Calculate transfer time
        function transfer_time = calculateTransferTime(obj)
            transfer_time = pi * sqrt(obj.a_transfer^3 / obj.MU);
        end
        
        % Print transfer details
        function printTransferDetails(obj)
            fprintf('\nHohmann Transfer Orbit Details:\n');
            fprintf('Initial Orbit Altitude: %.2f km\n', obj.r1 - obj.R_EARTH);
            fprintf('Final Orbit Altitude: %.2f km\n', obj.r2 - obj.R_EARTH);
            fprintf('Transfer Orbit Semi-Major Axis: %.2f km\n', obj.a_transfer);
            fprintf('Delta-V (Departure Burn): %.2f km/s\n', obj.delta_v_departure);
            fprintf('Delta-V (Arrival Burn): %.2f km/s\n', obj.delta_v_arrival);
            fprintf('Total Delta-V: %.2f km/s\n', obj.delta_v_departure + obj.delta_v_arrival);
            fprintf('Transfer Time: %.2f hours\n', obj.calculateTransferTime() / 3600);
        end
        
        % Visualize transfer orbit (similar to Python implementation)
        function visualizeTransfer(obj)
            % Create figure
            figure('Position', [100, 100, 1200, 800]);
            
            % 3D plot setup
            hold on;
            axis equal;
            
            % Earth representation (wireframe sphere)
            [x, y, z] = sphere(20);
            x = x * obj.R_EARTH;
            y = y * obj.R_EARTH;
            z = z * obj.R_EARTH;
            surf(x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'blue');
            
            % Orbit visualization
            theta = linspace(0, 2*pi, 100);
            
            % Initial orbit
            orbit1_x = obj.r1 * cos(theta);
            orbit1_y = obj.r1 * sin(theta);
            plot3(orbit1_x, orbit1_y, zeros(size(theta)), 'g-', 'LineWidth', 2, 'DisplayName', 'Initial Orbit');
            
            % Final orbit
            orbit2_x = obj.r2 * cos(theta);
            orbit2_y = obj.r2 * sin(theta);
            plot3(orbit2_x, orbit2_y, zeros(size(theta)), 'r-', 'LineWidth', 2, 'DisplayName', 'Final Orbit');
            
            % Labeling
            xlabel('X (km)');
            ylabel('Y (km)');
            zlabel('Z (km)');
            title('Hohmann Transfer Orbit');
            legend('show');
            
            hold off;
            view(45, 20);
        end
    end
    
    % Static method for user interaction
    methods (Static)
        function main()
            try
                % User inputs
                initial_alt = input('Enter initial orbit altitude (km): ');
                final_alt = input('Enter final orbit altitude (km): ');
                
                % Create Hohmann transfer object
                transfer = HohmannTransfer(initial_alt, final_alt);
                
                % Print transfer details
                transfer.printTransferDetails();
                
                % Visualize transfer
                transfer.visualizeTransfer();
            catch ME
                % Error handling
                fprintf('Error: %s\n', ME.message);
                fprintf('Please enter valid numerical inputs.\n');
            end
        end
    end
end

% Optional: Uncomment the line below to run directly
% HohmannTransfer.main();
