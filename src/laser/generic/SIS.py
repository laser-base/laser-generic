"""
Export required components for an SIS model.

Agents transition from Susceptible to Infectious upon infection and are infectious for a duration.
Agents transition from Infectious back to Susceptible upon recovery.
"""

from laser.generic.components import InfectiousIS as Infectious
from laser.generic.components import Susceptible
from laser.generic.components import TransmissionSI as Transmission
from laser.generic.shared import State

__all__ = ["Infectious", "State", "Susceptible", "Transmission"]
