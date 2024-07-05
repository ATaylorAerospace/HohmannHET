% Hohmann Orbit Transfer Analysis

% Constants
mu = 398600.4418; % Earth's gravitational parameter (km^3/s^2)
req = 6378.137; % Earth's equatorial radius (km)

% Conversion factor
rtd = 180/pi; % Radians to degrees
dtr = pi/180; % Degrees to radians

% Brent root-finding tolerance
rtol = 1e-8;

% Requesting user inputs and converting degrees to radians for orbit calculations
alt1 = input('Please input the initial altitude (kilometers): ');
alt2 = input('Please input the final altitude (kilometers): ');
inc1 = input('Please input the initial orbital inclination (degrees): ');
inc2 = input('Please input the final orbital inclination (degrees): ');
inc1 = inc1 * dtr;
inc2 = inc2 * dtr;

% Calculations for Hohmann transfer
dinc = abs(inc2 - inc1);
r1 = req + alt1;
r2 = req + alt2;
v1 = sqrt(mu / r1);
hn1 = sqrt(2 * r2 / (r1 + r2));
hn2 = sqrt(r1 / r2);
hn3 = sqrt(2 * r1 / (r1 + r2));

% Earth representation (3D plot)
[X, Y] = sphere(100);
R = req * ones(size(X));
fig = figure;
plot3(X, Y, R, 'Color', 'lightblue', 'Alpha', 0.6);
title('Hohmann Transfer: Earth and Orbits');
xlabel('X'); ylabel('Y'); zlabel('Z');

% Orbits visualization (simplified for demonstration)
theta = 0:0.01:2*pi; % Small angle increments for plotting
orbit1 = [r1*cos(theta), r1*sin(theta), zeros(size(theta))];
orbit2 = [r2*cos(theta), r2*sin(theta), zeros(size(theta))];

% Plot initial and final orbits
plot3(orbit1(1), orbit1(2), orbit1(3), 'LineWidth', 2);
hold on
plot3(orbit2(1), orbit2(2), orbit2(3), 'LineWidth', 2);

% Plot inclination change
dinc1 = inc1 - pi/2;
dinc2 = inc2 - pi/2;
plot3([0, r1*cos(dinc1)], [0, r1*sin(dinc1)], ...
    [0, hn1*r1*sin(dinc1/2)], 'LineWidth', 2, 'Color', 'red');
plot3([0, r2*cos(dinc2)], [0, r2*sin(dinc2)], ...
    [0, hn1*r2*sin((pi-dinc2)/2)], 'LineWidth', 2, 'Color', 'red');

% Labeling and grid
grid on;
view(-25, 30); % Adjusting view angle
legend('Initial Orbit', 'Final Orbit', 'Inclination Change');

% Display figure
show function delta_inc = hohm_inclination(alt1, alt2, inc1, inc2)
    % Calculate the change in inclination angle for Hohmann transfer
    % Inputs:
    %   alt1, alt2 (km): Initial and final altitudes
    %   inc1, inc2 (deg): Initial and final orbital inclinations
    % Output:
    %   delta_inc (rad): Change in inclination angle

    r1 = req + alt1;
    r2 = req + alt2;
    dinc = abs(inc2 - inc1);
    
    hn1 = sqrt(2 * r2 / (r1 + r2));
    hn2 = sqrt(r1 / r2);
    
    delta_inc = pi - 2*atan(hn1*sin(dinc/2));
end

% Find the delta inclination angle
delta_inc = hohm_inclination(alt1, alt2, inc1, inc2);
delta_inc_deg = delta_inc / dtr;

% Calculate the minimum delta v required for transfer
v1 = sqrt(mu / r1);
v2 = sqrt(mu / r2);
Delta_v = 2 * v1 * sin(delta_inc / 2);

% Display results
fprintf('Change in Inclination: %.2f degrees\n', delta_inc_deg);
fprintf('Minimum Delta-v: %.2f m/s\n', Delta_v);
```
