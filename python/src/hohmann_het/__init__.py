# Author: A Taylor | Purpose: Package initialisation | Ref: Vallado/Curtis
"""hohmann_het - Low-thrust orbital transfer library.

Modules
-------
dynamics     : Keplerian Hohmann mechanics, impulsive maneuvers, TOF.
propulsion   : High-fidelity Hall Effect Thruster (HET) models.
optimization : Solvers for minimum-propellant / minimum-TOF missions.
"""
from .dynamics import compute_hohmann, circular_velocity, orbital_period
from .propulsion import compute_het_state, propellant_mass, burn_time
from .optimization import optimize_isp, min_propellant_transfer

__version__ = "1.0.0"
__author__ = "A Taylor"
__all__ = [
    "compute_hohmann",
    "circular_velocity",
    "orbital_period",
    "compute_het_state",
    "propellant_mass",
    "burn_time",
    "optimize_isp",
    "min_propellant_transfer",
]
