import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi, g

# Calculation Constants
mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)
req = 6378.137  # Earth's equatorial radius (km)
rho_0 = 1.225  # Sea level air density (kg/m³)
H = 7.2  # Scale height of atmosphere (km)

# Conversion factors
rtd = 180 / pi  # Radians to degrees
dtr = pi / 180  # Degrees to radians

def atmospheric_density(altitude):
    """
    Calculate atmospheric density using exponential atmospheric model
    
    Args:
        altitude (float): Altitude above Earth's surface (km)
    
    Returns:
        float: Atmospheric density (kg/m³)
    """
    return rho_0 * np.exp(-altitude / H)

def drag_force(altitude, velocity, area, drag_coefficient=2.2):
    """
    Calculate atmospheric drag force
    
    Args:
        altitude (float): Altitude above Earth's surface (km)
        velocity (float): Orbital velocity (km/s)
        area (float): Cross-sectional area (m²)
        drag_coefficient (float): Drag coefficient (typically 2.2 for spacecraft)
    
    Returns:
        float: Drag force (N)
    """
    # Convert velocity to m/s
    v_ms = velocity * 1000
    
    # Calculate atmospheric density
    rho = atmospheric_density(altitude)
    
    # Drag force calculation
    drag = 0.5 * rho * (v_ms ** 2) * area * drag_coefficient
    
    return drag

def orbital_lifetime(altitude, velocity, mass, area, drag_coefficient=2.2):
    """
    Estimate orbital lifetime due to atmospheric drag
    
    Args:
        altitude (float): Initial altitude (km)
        velocity (float): Orbital velocity (km/s)
        mass (float): Spacecraft mass (kg)
        area (float): Cross-sectional area (m²)
        drag_coefficient (float): Drag coefficient
    
    Returns:
        float: Estimated orbital lifetime (days)
    """
    # Gravitational and orbital parameters
    r = req + altitude
    
    # Orbital period
    period = 2 * pi * np.sqrt(r**3 / mu)
    
    # Drag force
    drag = drag_force(altitude, velocity, area, drag_coefficient)
    
    # Power loss calculation
    power_loss = drag * velocity * 1000  # Convert velocity to m/s
    
    # Calculate orbital decay rate
    orbital_decay_rate = 2 * power_loss / (mass * g * period * 1000)  # Convert period to seconds
    
    # Estimate lifetime before re-entry
    altitude_loss_per_orbit = orbital_decay_rate * period * req
    total_lifetime = altitude / altitude_loss_per_orbit * period / (24 * 3600)  # Convert to days
    
    return total_lifetime

def calculate_orbital_parameters(altitude):
    """
    Calculate basic orbital parameters
    
    Args:
        altitude (float): Altitude above Earth's surface (km)
    
    Returns:
        tuple: Orbital radius, orbital velocity, orbital period
    """
    r = req + altitude
    
    # Orbital velocity
    v = np.sqrt(mu / r)
    
    # Orbital period
    period = 2 * pi * np.sqrt(r**3 / mu)
    
    return r, v, period

def main():
    # User inputs for spacecraft parameters
    print("\nOrbital Drag Analysis")
    
    # Orbital parameters
    alt = float(input("Enter orbital altitude (km): "))
    mass = float(input("Enter spacecraft mass (kg): "))
    area = float(input("Enter cross-sectional area (m²): "))
    
    # Calculate orbital parameters
    r, v, period = calculate_orbital_parameters(alt)
    
    # Atmospheric density and drag calculations
    rho = atmospheric_density(alt)
    drag = drag_force(alt, v, area)
    lifetime = orbital_lifetime(alt, v, mass, area)
    
    # Output results
    print("\nOrbital Parameters:")
    print(f"Orbital Radius: {r:.2f} km")
    print(f"Orbital Velocity: {v:.2f} km/s")
    print(f"Orbital Period: {period/60:.2f} minutes")
    print("\nAtmospheric Drag Analysis:")
    print(f"Atmospheric Density: {rho:.6f} kg/m³")
    print(f"Drag Force: {drag:.4f} N")
    print(f"Estimated Orbital Lifetime: {lifetime:.2f} days")
    
    # Visualization of drag effects
    altitudes = np.linspace(200, 400, 100)
    densities = [atmospheric_density(h) for h in altitudes]
    
    plt.figure(figsize=(10, 6))
    plt.plot(altitudes, densities)
    plt.title('Atmospheric Density vs Altitude')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Atmospheric Density (kg/m³)')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
