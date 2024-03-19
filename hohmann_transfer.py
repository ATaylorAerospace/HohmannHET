import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi
from scipy.optimize import brentq
from mpl_toolkits.mplot3d import Axes3D

# Constants
mu = 398600.4418  # km^3/s^2, Earth's gravitational parameter
req = 6378.137  # km, Earth's equatorial radius
rtd = 180 / np.pi
dtr = np.pi / 180

# Brent root-finding tolerance
rtol = 1.0e-8

def get_input(prompt, condition=lambda x: True):
    while True:
        try:
            value = float(input(prompt))
            if condition(value):
                return value
        except ValueError:
            pass

def hohmfunc(dinc1, v1, hn1, hn2, hn3, dinc):
    dinc2 = dinc - dinc1
    return v1 * np.sqrt(1.0 + hn1**2 - 2 * hn1 * np.cos(dinc1)) + \
           v1 * np.sqrt(hn2**2 * hn3**2 + hn2**2 - 2 * hn2**2 * hn3 * np.cos(dinc2))

# Requesting user inputs
print('\nHohmann Orbit Transfer Analysis\n')
alt1 = get_input('Please input the initial altitude (kilometers): ', lambda x: x > 0)
alt2 = get_input('Please input the final altitude (kilometers): ', lambda x: x > 0)
inc1 = get_input('Please input the initial orbital inclination (degrees) (0 <= inclination <= 180): ', lambda x: 0 <= x <= 180) * dtr
inc2 = get_input('Please input the final orbital inclination (degrees) (0 <= inclination <= 180): ', lambda x: 0 <= x <= 180) * dtr

# Calculating transfer parameters
dinc = abs(inc2 - inc1)
r1, r2 = req + alt1, req + alt2
hn1, hn2, hn3 = np.sqrt(2 * r2 / (r1 + r2)), np.sqrt(r1 / r2), np.sqrt(2 * r1 / (r1 + r2))
v1, v2 = np.sqrt(mu / r1), np.sqrt(mu / r2)
smat = 0.5 * (r1 + r2)
ecct = (max(r1, r2) - min(r1, r2)) / (r1 + r2)
rp, ra = smat * (1 - ecct), smat * (1 + ecct)
vt1, vt2 = np.sqrt(2 * mu * ra / (rp * (ra + rp))), np.sqrt(2 * mu * rp / (ra * (ra + rp)))
taut, tof = 2 * pi * np.sqrt(smat**3 / mu), pi * np.sqrt(smat**3 / mu)

# Non-coplanar transfer adjustments
if dinc == 0:
    dv1, dv2 = np.abs(vt1 - v1), np.abs(vt2 - v2)
    dinc1, dinc2, inct = 0, 0, inc1
else:
    # Root finding for dinc1, then calculate dinc2
    xroot = brentq(hohmfunc, 0, dinc, args=(v1, hn1, hn2, hn3, dinc), xtol=rtol)
    dinc1, dinc2 = xroot, dinc - xroot
    dv1 = v1 * np.sqrt(1.0 + hn1**2 - 2 * hn1 * np.cos(dinc1))
    dv2 = v1 * np.sqrt(hn2**2 * hn3**2 + hn2**2 - 2 * hn2**2 * hn3 * np.cos(dinc2))
    inct = inc1 + dinc1 if inc2 > inc1 else inc1 - dinc1

# Output results
print(f'\nHohmann Orbit Transfer Analysis\n{"-"*30}\n')
print(f'Initial orbit altitude: {alt1:.4f} km')
# Additional print statements for each calculated parameter similar to the MATLAB script

# Plotting orbits and primer vector analysis
# This would typically require converting state vectors to position and plotting them
# For simplicity, we can plot spheres for Earth and circular orbits for illustration

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Earth
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
earth_x = earth_radius * np.cos(u) * np.sin(v)
earth_y = earth_radius * np.sin(u) * np.sin(v)
earth_z = earth_radius * np.cos(v)
ax.plot_surface(earth_x, earth_y, earth_z, color='lightblue', alpha=0.6)

# Placeholder orbits (example, not accurate)
theta = np.linspace(0, 2 * np.pi, 100)
orbit1_x = (r1 * np.cos(theta)) / req
orbit1_y = (r1 * np.sin(theta)) / req
orbit2_x = (r2 * np.cos(theta)) / req
orbit2_y = (r2 * np.sin(theta)) / req

ax.plot(orbit1_x, orbit1_y, zs=0, zdir='z', label='Initial Orbit')
ax.plot(orbit2_x, orbit2_y, zs=0, zdir='z', label='Final Orbit')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('Hohmann Transfer: Initial and Final Orbits')
plt.legend()
plt.show()
