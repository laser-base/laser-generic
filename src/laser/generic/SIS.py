"""
Export required components for an SIS model.

Agents transition from Susceptible to Infectious upon infection and are infectious for a duration.
Agents transition from Infectious back to Susceptible upon recovery.
"""

from .components import InfectiousIS
from .components import Susceptible
from .components import TransmissionSI
from .shared import State

__all__ = ["InfectiousIS", "State", "Susceptible", "TransmissionSI"]
