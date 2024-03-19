# Hohmann Orbit Transfer Analysis

This script performs a Hohmann transfer orbit analysis, including calculations for both planar and non-coplanar circular orbit transfers. It includes generating three-dimensional orbit graphics and a graphical primer vector analysis.

## Overview

The Hohmann transfer is the most efficient way to move a spaceprobes between two orbits using two engine impulses. This script provides a tool to analyze such transfers, calculating key parameters like delta-v (change in velocity), transfer orbit characteristics, and inclination changes.

## Features

- Calculation of transfer parameters for orbits around Earth.
- Support for planar and non-coplanar transfers.
- 3D visualization of the initial, transfer, and final orbits.
- Interactive input for defining initial and final orbit parameters.

## Dependencies

- `numpy`
- `matplotlib`
- `scipy`

Ensure these packages are installed in your Python environment to run the script successfully.

## Usage

1. Clone this repository or download the script.
2. Ensure you have Python 3.x installed along with the required dependencies.
3. Run the script using a Python interpreter:
   ```
   python hohmann_transfer.py
   ```
4. Follow the on-screen prompts to input the initial altitude, final altitude, initial orbital inclination, and final orbital inclination.
5. The script will calculate and display the transfer analysis results, followed by a 3D visualization of the orbits.

## Mathematical Background

The script employs several key concepts from celestial mechanics:

- **Gravitational Parameter (\(μ\))**: Earth's gravitational constant times the mass of Earth.
- **Radius of Earth (\(R_{eq}\))**: Equatorial radius of Earth.
- **Delta-v (\(\Delta v\))**: The measure of impulse required to perform the transfer.
- **Semi-major axis of the transfer orbit**: Calculated as the average of the initial and final orbit radii.
- **Eccentricity of the transfer orbit**: Determined based on the radii of the initial and final orbits.

## Visualization

The 3D visualization plots the Earth, initial orbit, and final orbit to provide a graphical understanding of the transfer. Note: The visualization is simplified for educational purposes and might not represent exact scales.

## Contributing


## License

This project is open-sourced under the MIT License. See the [LICENSE](LICENSE) file for more details.

---
