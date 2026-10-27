# Author: A Taylor | Purpose: Package initialisation | Ref: Vallado/Curtis
"""hohmann_het - Low-thrust orbital transfer library.

Modules
-------
dynamics     : Keplerian Hohmann mechanics, impulsive maneuvers, TOF.
propulsion   : High-fidelity Hall Effect Thruster (HET) models.
optimization : Solvers for minimum-propellant / minimum-TOF missions.
"""
from .dynamics import (
    MU,
    R_EARTH,
    HohmannTransfer,
    circular_velocity,
    compute_hohmann,
    orbital_period,
)
from .optimization import (
    OptimizationResult,
    golden_section_minimize,
    min_propellant_transfer,
    optimize_isp,
)
from .propulsion import (
    E_CHARGE,
    G0,
    M_XE,
    HETState,
    burn_time,
    compute_het_state,
    propellant_mass,
)

__version__ = "1.0.0"
__author__ = "A Taylor"
__all__ = [
    "E_CHARGE",
    "G0",
    "MU",
    "M_XE",
    "R_EARTH",
    "HETState",
    "HohmannTransfer",
    "OptimizationResult",
    "burn_time",
    "circular_velocity",
    "compute_het_state",
    "compute_hohmann",
    "golden_section_minimize",
    "min_propellant_transfer",
    "optimize_isp",
    "orbital_period",
    "propellant_mass",
]
