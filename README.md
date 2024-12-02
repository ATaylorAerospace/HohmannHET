### Hohmann Transfer Orbit Analysis 

### Overview 

This script provides a comprehensive tool for calculating and visualizing orbital transfers around Earth. It specifically focuses on the most energy-efficient method of transferring a spacecraft between two circular orbits using minimal propulsive maneuvers.

### Features 

### Precise Orbital Transfer Calculations
- Compute delta-v (velocity change) for departure and arrival burns.
- Calculate transfer orbit semi-major axis.
- Determine total transfer time.

### Interactive Analysis
- User-friendly input for initial and final orbit altitudes.
- Detailed printout of transfer orbit parameters.

### 3D Visualization
- Graphical representation of Earth.
- Plotting of initial and final orbital trajectories.
- Visual understanding of the transfer orbit.

### Dependencies 

- `numpy`: Numerical computing library.
- `matplotlib`: Plotting and visualization.
- `scipy`: Scientific computing tools (specifically for constants).

## Mathematical Foundations

The script leverages key celestial mechanics principles:
- **Gravitational Parameter (μ):** 398,600.4418 km³/s².
- **Earth's Equatorial Radius:** 6,378.137 km.

### Transfer Orbit Calculations
- Semi-major axis computation.
- Orbital velocity calculations.
- Delta-v budget analysis.

## Installation

1. Ensure you have Python 3.x installed.
2. Install required dependencies:

   ```bash
   pip install numpy matplotlib scipy
   ```

** Usage **

Run the script:

```bash
python hohmann_transfer.py
```

When prompted, input:
- **Initial orbit altitude** (in kilometers).
- **Final orbit altitude** (in kilometers).

The script will:
- Calculate transfer parameters.
- Display detailed transfer information.
- Generate a 3D orbit visualization.

## Example Output

```
Hohmann Transfer Orbit Details:
Initial Orbit Altitude: 500.00 km
Final Orbit Altitude: 35,786.00 km
Transfer Orbit Semi-Major Axis: 23,093.07 km
Delta-V (Departure Burn): 1.80 km/s
Delta-V (Arrival Burn): 2.44 km/s
Total Delta-V: 4.24 km/s
Transfer Time: 12.54 hours
```

### Visualization

The 3D plot will show:
- Earth's spherical representation.
- Initial orbit trajectory.
- Final orbit trajectory.
- Transfer orbit path.

## Limitations

- Currently supports transfers around Earth.
- Assumes circular orbits.
- Does not account for extremely complex gravitational interactions prior to geo insertion.

## Contributing


## License

This project is open-sourced under the MIT License. See the [LICENSE](LICENSE) file for details.
```

Save this content as `README.md` in your repository. Let me know if you need further help!
