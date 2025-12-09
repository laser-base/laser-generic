"""
Export required components for an SIRS model.

Agents transition from Susceptible to Infectious upon infection.
Agents transition from Infectious to Recovered upon recovery after the infectious duration.
Agents transition from Recovered back to Susceptible upon waning immunity after the waning duration.
"""

from .components import InfectiousIRS
from .components import RecoveredRS
from .components import Susceptible
from .components import TransmissionSI
from .shared import State

__all__ = ["InfectiousIRS", "RecoveredRS", "State", "Susceptible", "TransmissionSI"]
