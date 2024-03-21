import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi
from scipy.optimize import brentq

# Constants
mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)
req = 6378.137  # Earth's equatorial radius (km)

# Conversion factors
rtd = 180 / pi  # Radians to degrees
dtr = pi / 180  # Degrees to radians

# Brent root-finding tolerance
rtol = 1.0e-8

def get_input(prompt, condition=lambda x: True):
    while True:
        try:
            value = float(input(prompt))
            if condition(value):
                return value
        except ValueError:
            print("Invalid input. Please try again.")
            continue

def hohmfunc(dinc1, v1, hn1, hn2, hn3, dinc):
    dinc2 = dinc - dinc1
    return v1 * np.sqrt(1.0 + hn1**2 - 2 * hn1 * np.cos(dinc1 * dtr)) + \
           v1 * np.sqrt(hn2**2 * hn3**2 + hn2**2 - 2 * hn2**2 * hn3 * np.cos(dinc2 * dtr))

# Requesting user inputs and converting degrees to radians for calculations
print('\nHohmann Orbit Transfer Analysis\n')
alt1 = get_input('Please input the initial altitude (kilometers): ', lambda x: x > 0)
alt2 = get_input('Please input the final altitude (kilometers): ', lambda x: x > 0)
inc1 = get_input('Please input the initial orbital inclination (degrees) (0 <= inclination <= 180): ', lambda x: 0 <= x <= 180) * dtr
inc2 = get_input('Please input the final orbital inclination (degrees) (0 <= inclination <= 180): ', lambda x: 0 <= x <= 180) * dtr

# Calculations for Hohmann transfer
dinc = abs(inc2 - inc1)
r1, r2 = req + alt1, req + alt2
v1 = np.sqrt(mu / r1)
hn1, hn2, hn3 = np.sqrt(2 * r2 / (r1 + r2)), np.sqrt(r1 / r2), np.sqrt(2 * r1 / (r1 + r2))

# Plotting preparations
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Earth representation
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = req * np.cos(u) * np.sin(v)
y = req * np.sin(u) * np.sin(v)
z = req * np.cos(v)
ax.plot_surface(x, y, z, color='lightblue', alpha=0.6)

# Orbits visualization (simplified for demonstration)
theta = np.linspace(0, 2 * np.pi, 100)
orbit1_x, orbit1_y = (r1 * np.cos(theta)), (r1 * np.sin(theta))
orbit2_x, orbit2_y = (r2 * np.cos(theta)), (r2 * np.sin(theta))

ax.plot(orbit1_x, orbit1_y, zs=0, zdir='z', label='Initial Orbit')
ax.plot(orbit2_x, orbit2_y, zs=0, zdir='z', label='Final Orbit')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('Hohmann Transfer: Initial and Final Orbits')
plt.legend()
plt.show()
