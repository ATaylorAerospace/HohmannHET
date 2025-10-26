###  🚀 Hohmann Transfer with Hall Effect Thrusters (HohmannHET)

[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Stars](https://img.shields.io/github/stars/ATaylorAerospace/HohmannHET?style=social)](https://github.com/ATaylorAerospace/HohmannHET)
[![Language](https://img.shields.io/badge/Languages-Python%20%7C%20C%2B%2B%20%7C%20MATLAB-brightgreen)](https://github.com/ATaylorAerospace/HohmannHET)

## 📋 Overview

**HohmannHET** is a comprehensive multi-language toolkit for calculating and visualizing Hohmann transfer orbits around Earth. This repository provides highly optimized implementations in **Python**, **C++**, and **MATLAB**, specifically designed for analyzing the most energy-efficient orbital transfer maneuvers between circular orbits.

The project focuses on precision, performance, and educational value, making it suitable for both academic research and practical aerospace applications.

---

## ✨ Key Features

### 🎯 **Precise Orbital Mechanics Calculations**
- **Delta-V Computation**: Calculate accurate velocity changes for departure and arrival burns
- **Transfer Orbit Analysis**: Determine semi-major axis and orbital parameters
- **Time Optimization**: Compute total transfer duration with high precision
- **Multi-Scenario Support**: Handle LEO-to-GEO, reverse transfers, and custom altitudes

### 🔧 **Multi-Language Implementation**
- **Python**: Ultra-optimized with NumPy vectorization and matplotlib visualization
- **C++**: High-performance implementation with comprehensive Google Test suite
- **MATLAB**: Educational-friendly version for academic environments

### 📊 **Advanced Visualization**
- **3D Orbital Plots**: Interactive Earth representation with orbital trajectories
- **Transfer Path Visualization**: Clear depiction of elliptical transfer orbits
- **Performance Optimized**: Reduced rendering overhead for faster plotting

### 🧪 **Testing & Validation**
- **Comprehensive Test Suite**: 15+ test cases covering edge cases and boundary conditions
- **Performance Benchmarks**: Sub-millisecond calculation times
- **Mathematical Validation**: Verified against known LEO-to-GEO transfer values

---

## 🛠️ Technical Specifications

### **Physical Constants**
| Parameter | Value | Unit |
|-----------|-------|------|
| Earth's Gravitational Parameter (μ) | 398,600.4418 | km³/s² |
| Earth's Equatorial Radius | 6,378.137 | km |
| Calculation Precision | 1e-6 | - |

### **Supported Transfer Types**
- ✅ Low Earth Orbit (LEO) to Geostationary Orbit (GEO)
- ✅ Reverse transfers (GEO to LEO)
- ✅ Custom altitude transfers
- ✅ Same-altitude validation (zero delta-V)
- ✅ Extreme altitude differences (up to 100,000 km)

---

## 🚀 Quick Start

#### Prerequisites

```bash
pip install numpy matplotlib scipy
```

#### Citations
```
@misc{ATaylor_HohmannTransferOrbitAnalysis_2024,
  author       = {A. Taylor},
  title        = {Hohmann Transfer Orbit Analysis},
  year         = {2025},
  url          = {[https://github.com/ATaylorAerospace/HohmannHET]},
  note         = {Accessed: YYYY-MM-DD}
}
```
