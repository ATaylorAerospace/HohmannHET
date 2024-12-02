import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi

# Calculation Constants
MU = 398600.4418  # Earth's gravitational parameter (km^3/s^2)
R_EARTH = 6378.137  # Earth's equatorial radius (km)

class HohmannTransfer:
    def __init__(self, initial_altitude, final_altitude):
        """
        Initialize Hohmann Transfer Orbit Calculation
        
        Args:
            initial_altitude (float): Initial orbit altitude (km)
            final_altitude (float): Final orbit altitude (km)
        """
        # Validate inputs
        if initial_altitude <= 0 or final_altitude <= 0:
            raise ValueError("Altitudes must be positive")
        
        # Orbital radii
        self.r1 = R_EARTH + initial_altitude
        self.r2 = R_EARTH + final_altitude
        
        # Calculate transfer orbit parameters
        self.a_transfer = (self.r1 + self.r2) / 2  # Transfer orbit semi-major axis
        
        # Initial orbital velocity
        self.v1 = np.sqrt(MU / self.r1)
        
        # Compute delta-V for transfer
        self.delta_v_departure = self.calculate_delta_v_departure()
        self.delta_v_arrival = self.calculate_delta_v_arrival()
    
    def calculate_delta_v_departure(self):
        """
        Calculate delta-V for departure burn
        
        Returns:
            float: Delta-V for departure (km/s)
        """
        v_transfer_1 = np.sqrt(MU * (2/self.r1 - 1/self.a_transfer))
        return v_transfer_1 - self.v1
    
    def calculate_delta_v_arrival(self):
        """
        Calculate delta-V for arrival burn
        
        Returns:
            float: Delta-V for arrival (km/s)
        """
        # Velocity at final orbit
        v2 = np.sqrt(MU / self.r2)
        
        # Velocity at transfer orbit apoapsis
        v_transfer_2 = np.sqrt(MU * (2/self.r2 - 1/self.a_transfer))
        
        return v2 - v_transfer_2
    
    def calculate_transfer_time(self):
        """
        Calculate Hohmann transfer orbit time
        
        Returns:
            float: Transfer time (seconds)
        """
        return np.pi * np.sqrt(self.a_transfer**3 / MU)
    
    def visualize_transfer(self):
        """
        Visualize Hohmann transfer orbit
        """
        # Create figure
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Earth representation
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = R_EARTH * np.cos(u) * np.sin(v)
        y = R_EARTH * np.sin(u) * np.sin(v)
        z = R_EARTH * np.cos(v)
        ax.plot_surface(x, y, z, color='blue', alpha=0.2)
        
        # Orbit visualization
        theta = np.linspace(0, 2*np.pi, 100)
        
        # Initial orbit
        orbit1_x = self.r1 * np.cos(theta)
        orbit1_y = self.r1 * np.sin(theta)
        ax.plot(orbit1_x, orbit1_y, zs=0, zdir='z', label='Initial Orbit', color='green')
        
        # Final orbit
        orbit2_x = self.r2 * np.cos(theta)
        orbit2_y = self.r2 * np.sin(theta)
        ax.plot(orbit2_x, orbit2_y, zs=0, zdir='z', label='Final Orbit', color='red')
        
        # Labeling
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
        plt.title('Hohmann Transfer Orbit')
        plt.legend()
        plt.show()
    
    def print_transfer_details(self):
        """
        Print detailed information about the Hohmann transfer
        """
        print("\nHohmann Transfer Orbit Details:")
        print(f"Initial Orbit Altitude: {self.r1 - R_EARTH:.2f} km")
        print(f"Final Orbit Altitude: {self.r2 - R_EARTH:.2f} km")
        print(f"Transfer Orbit Semi-Major Axis: {self.a_transfer:.2f} km")
        print(f"Delta-V (Departure Burn): {self.delta_v_departure:.2f} km/s")
        print(f"Delta-V (Arrival Burn): {self.delta_v_arrival:.2f} km/s")
        print(f"Total Delta-V: {self.delta_v_departure + self.delta_v_arrival:.2f} km/s")
        print(f"Transfer Time: {self.calculate_transfer_time() / 3600:.2f} hours")

def main():
    """
    Main function to run Hohmann transfer analysis
    """
    try:
        # User inputs
        initial_alt = float(input("Enter initial orbit altitude (km): "))
        final_alt = float(input("Enter final orbit altitude (km): "))
        
        # Create Hohmann transfer object
        transfer = HohmannTransfer(initial_alt, final_alt)
        
        # Print transfer details
        transfer.print_transfer_details()
        
        # Visualize transfer
        transfer.visualize_transfer()
    
    except ValueError as e:
        print(f"Error: {e}")
        print("Please enter valid numerical inputs.")

if __name__ == "__main__":
    main()
