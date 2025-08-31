import numpy as np
from typing import Tuple, Optional
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

# Module-level constants for better performance (avoid dataclass overhead)
MU = 398600.4418      # Earth's gravitational parameter (km^3/s^2)
R_EARTH = 6378.137    # Earth's equatorial radius (km)
INV_3600 = 1.0 / 3600.0  # Precomputed hours conversion constant

# Precomputed visualization constants
SPHERE_PHI_POINTS = 20    # Reduced from 30 for faster rendering
SPHERE_THETA_POINTS = 15  # Reduced from 20 for faster rendering
ORBIT_POINTS = 80         # Reduced from 100 for faster plotting

class HohmannTransferOptimized:
    """Ultra-efficient implementation of Hohmann transfer calculations"""
    
    __slots__ = ('r1', 'r2', 'a_transfer', 'v1', 'delta_v_departure', 
                 'delta_v_arrival', 'total_delta_v', 'transfer_time', 
                 'transfer_hours', '_initial_alt', '_final_alt',
                 '_cos_theta', '_sin_theta', '_zeros')
    
    def __init__(self, initial_altitude: float, final_altitude: float):
        """
        Initialize Hohmann Transfer with optimized calculations
        
        Args:
            initial_altitude: Initial orbit altitude (km)
            final_altitude: Final orbit altitude (km)
        """
        # Fast validation without exception overhead for valid inputs
        if initial_altitude <= 0 or final_altitude <= 0:
            raise ValueError("Altitudes must be positive")
        
        # Store original altitudes for display
        self._initial_alt = initial_altitude
        self._final_alt = final_altitude
        
        # Precompute all orbital parameters
        self.r1 = R_EARTH + initial_altitude
        self.r2 = R_EARTH + final_altitude
        self.a_transfer = (self.r1 + self.r2) * 0.5
        
        # Vectorized velocity calculations
        self._compute_all_values()
        
        # Precompute visualization arrays
        self._precompute_visualization()

    def _compute_all_values(self) -> None:
        """Compute all orbital parameters in one vectorized pass"""
        # Vectorized calculations using numpy arrays
        radii = np.array([self.r1, self.r2])
        mu_over_r = MU / radii
        v_circular = np.sqrt(mu_over_r)  # [v1, v2]
        
        self.v1 = v_circular[0]
        
        # Transfer orbit velocities using optimized calculations
        mu_2_over_r = 2.0 * mu_over_r
        mu_over_a = MU / self.a_transfer
        
        v_transfer = np.sqrt(mu_2_over_r - mu_over_a)  # [v_transfer_1, v_transfer_2]
        
        # Delta-V calculations
        self.delta_v_departure = v_transfer[0] - v_circular[0]
        self.delta_v_arrival = v_circular[1] - v_transfer[1]
        self.total_delta_v = self.delta_v_departure + self.delta_v_arrival
        
        # Precompute transfer time
        self.transfer_time = np.pi * np.sqrt(self.a_transfer**3 / MU)
        self.transfer_hours = self.transfer_time * INV_3600

    def _precompute_visualization(self) -> None:
        """Precompute trigonometric values for visualization"""
        theta = np.linspace(0, 2*np.pi, ORBIT_POINTS)
        self._cos_theta = np.cos(theta)
        self._sin_theta = np.sin(theta)
        self._zeros = np.zeros_like(theta)

    def visualize_transfer(self) -> Figure:
        """
        Visualize Hohmann transfer orbit with maximum optimization
        
        Returns:
            matplotlib.figure.Figure: The generated figure
        """
        # Create figure with optimized backend
        fig = plt.figure(figsize=(12, 8))
        ax: Axes = fig.add_subplot(111, projection='3d')
        
        # Ultra-optimized Earth sphere using meshgrid
        phi = np.linspace(0, 2*np.pi, SPHERE_PHI_POINTS)
        theta = np.linspace(0, np.pi, SPHERE_THETA_POINTS)
        phi_mesh, theta_mesh = np.meshgrid(phi, theta)
        
        # Vectorized sphere coordinates
        sin_theta = np.sin(theta_mesh)
        x = R_EARTH * np.cos(phi_mesh) * sin_theta
        y = R_EARTH * np.sin(phi_mesh) * sin_theta
        z = R_EARTH * np.cos(theta_mesh)
        
        # Plot Earth with minimal overhead
        ax.plot_surface(x, y, z, color='blue', alpha=0.15, 
                       linewidth=0, antialiased=False, shade=False)
        
        # Plot orbits using precomputed trigonometric values
        ax.plot(self.r1 * self._cos_theta, self.r1 * self._sin_theta, self._zeros, 
               label='Initial Orbit', color='green', linewidth=2)
        ax.plot(self.r2 * self._cos_theta, self.r2 * self._sin_theta, self._zeros,
               label='Final Orbit', color='red', linewidth=2)
        
        # Add transfer orbit ellipse
        self._plot_transfer_orbit(ax)
        
        # Optimized axis settings
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
        ax.set_title('Hohmann Transfer Orbit')
        ax.legend()
        ax.view_init(30, 45)
        
        return fig

    def _plot_transfer_orbit(self, ax: Axes) -> None:
        """Plot transfer orbit ellipse efficiently"""
        # Semi-ellipse for transfer orbit
        theta_transfer = np.linspace(0, np.pi, ORBIT_POINTS // 2)
        
        # Ellipse parameters
        a = self.a_transfer
        c = abs(self.r2 - self.r1) * 0.5
        b = np.sqrt(a*a - c*c) if a > c else a * 0.1  # Avoid complex numbers
        
        # Parametric ellipse
        x_ellipse = a * np.cos(theta_transfer)
        y_ellipse = b * np.sin(theta_transfer)
        
        # Center and plot
        x_center = (self.r1 + self.r2) * 0.5
        x_transfer = x_ellipse + x_center - a
        
        ax.plot(x_transfer, y_ellipse, np.zeros_like(theta_transfer),
               '--', color='magenta', linewidth=1.5, label='Transfer Orbit')

    def print_transfer_details(self) -> None:
        """Print transfer details with optimized string formatting"""
        # Single f-string for maximum efficiency
        print(f"""
Hohmann Transfer Orbit Details:
Initial Orbit Altitude: {self._initial_alt:.2f} km
Final Orbit Altitude: {self._final_alt:.2f} km
Transfer Orbit Semi-Major Axis: {self.a_transfer:.2f} km
Delta-V (Departure Burn): {self.delta_v_departure:.2f} km/s
Delta-V (Arrival Burn): {self.delta_v_arrival:.2f} km/s
Total Delta-V: {self.total_delta_v:.2f} km/s
Transfer Time: {self.transfer_hours:.2f} hours""")

    # Property-style getters for better performance than methods
    @property
    def departure_delta_v(self) -> float:
        return self.delta_v_departure
    
    @property
    def arrival_delta_v(self) -> float:
        return self.delta_v_arrival
    
    @property
    def total_delta_v_property(self) -> float:
        return self.total_delta_v
    
    @property
    def transfer_time_hours(self) -> float:
        return self.transfer_hours

def get_valid_input(prompt: str) -> float:
    """Optimized input validation with minimal overhead"""
    while True:
        try:
            value = float(input(prompt))
            if value > 0:
                return value
            print("Please enter a positive number.")
        except ValueError:
            print("Please enter a valid number.")

def main() -> None:
    """Optimized main function"""
    try:
        # Fast input with optimized validation
        initial_alt = get_valid_input("Enter initial orbit altitude (km): ")
        final_alt = get_valid_input("Enter final orbit altitude (km): ")
        
        # Create transfer object
        transfer = HohmannTransferOptimized(initial_alt, final_alt)
        
        # Display results
        transfer.print_transfer_details()
        
        # Visualize (optional - comment out for pure calculation performance)
        fig = transfer.visualize_transfer()
        plt.show()
        
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()
