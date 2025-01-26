import numpy as np
from dataclasses import dataclass
from typing import Tuple
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

# Constants as a frozen dataclass for better organization and performance
@dataclass(frozen=True)
class OrbitalConstants:
    MU: float = 398600.4418      # Earth's gravitational parameter (km^3/s^2)
    R_EARTH: float = 6378.137    # Earth's equatorial radius (km)

class HohmannTransfer:
    """Efficient implementation of Hohmann transfer calculations"""
    
    _constants = OrbitalConstants()
    
    def __init__(self, initial_altitude: float, final_altitude: float):
        """
        Initialize Hohmann Transfer Orbit Calculator
        
        Args:
            initial_altitude: Initial orbit altitude (km)
            final_altitude: Final orbit altitude (km)
        Raises:
            ValueError: If altitudes are not positive
        """
        if initial_altitude <= 0 or final_altitude <= 0:
            raise ValueError("Altitudes must be positive")
        
        # Precompute all orbital parameters at initialization
        self.r1 = self._constants.R_EARTH + initial_altitude
        self.r2 = self._constants.R_EARTH + final_altitude
        self.a_transfer = (self.r1 + self.r2) * 0.5
        
        # Precompute velocities and delta-Vs
        self.v1 = np.sqrt(self._constants.MU / self.r1)
        self.delta_v_departure, self.delta_v_arrival = self._calculate_delta_vs()
        self.total_delta_v = self.delta_v_departure + self.delta_v_arrival
        
        # Precompute transfer time
        self.transfer_time = np.pi * np.sqrt(self.a_transfer**3 / self._constants.MU)

    def _calculate_delta_vs(self) -> Tuple[float, float]:
        """Calculate both delta-Vs in one efficient computation"""
        # Use cached values and compute everything at once
        v_transfer_1 = np.sqrt(self._constants.MU * (2/self.r1 - 1/self.a_transfer))
        v_transfer_2 = np.sqrt(self._constants.MU * (2/self.r2 - 1/self.a_transfer))
        v2 = np.sqrt(self._constants.MU / self.r2)
        
        return v_transfer_1 - self.v1, v2 - v_transfer_2

    def visualize_transfer(self) -> Figure:
        """
        Visualize Hohmann transfer orbit with optimized plotting
        
        Returns:
            matplotlib.figure.Figure: The generated figure
        """
        fig = plt.figure(figsize=(12, 8))
        ax: Axes = fig.add_subplot(111, projection='3d')
        
        # Optimize Earth sphere representation
        phi = np.linspace(0, 2*np.pi, 30)
        theta = np.linspace(0, np.pi, 20)
        phi, theta = np.meshgrid(phi, theta)
        
        r = self._constants.R_EARTH
        x = r * np.cos(phi) * np.sin(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(theta)
        
        # Plot with optimized parameters
        ax.plot_surface(x, y, z, color='blue', alpha=0.2, shade=False)
        
        # Optimize orbit plotting using vectorization
        theta = np.linspace(0, 2*np.pi, 100)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        
        # Plot orbits efficiently
        zeros = np.zeros_like(theta)
        ax.plot(self.r1 * cos_theta, self.r1 * sin_theta, zeros, 
               label='Initial Orbit', color='green', linewidth=2)
        ax.plot(self.r2 * cos_theta, self.r2 * sin_theta, zeros,
               label='Final Orbit', color='red', linewidth=2)
        
        # Set visualization properties
        ax.set(xlabel='X (km)', ylabel='Y (km)', zlabel='Z (km)',
               title='Hohmann Transfer Orbit')
        ax.legend()
        ax.view_init(30, 45)
        
        return fig

    def print_transfer_details(self) -> None:
        """Print transfer details using efficient string formatting"""
        details = f"""
Hohmann Transfer Orbit Details:
Initial Orbit Altitude: {self.r1 - self._constants.R_EARTH:.2f} km
Final Orbit Altitude: {self.r2 - self._constants.R_EARTH:.2f} km
Transfer Orbit Semi-Major Axis: {self.a_transfer:.2f} km
Delta-V (Departure Burn): {self.delta_v_departure:.2f} km/s
Delta-V (Arrival Burn): {self.delta_v_arrival:.2f} km/s
Total Delta-V: {self.total_delta_v:.2f} km/s
Transfer Time: {self.transfer_time / 3600:.2f} hours"""
        print(details)

def main() -> None:
    """Main function with improved error handling"""
    try:
        initial_alt = float(input("Enter initial orbit altitude (km): "))
        final_alt = float(input("Enter final orbit altitude (km): "))
        
        transfer = HohmannTransfer(initial_alt, final_alt)
        transfer.print_transfer_details()
        
        fig = transfer.visualize_transfer()
        plt.show()
        
    except ValueError as e:
        print(f"Error: {str(e)}")
        print("Please enter valid numerical inputs.")
    except Exception as e:
        print(f"Unexpected error: {str(e)}")

if __name__ == "__main__":
    main()
