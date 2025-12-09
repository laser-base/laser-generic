"""
Export required components for an SIRS model.

Agents transition from Susceptible to Infectious upon infection.
Agents transition from Infectious to Recovered upon recovery after the infectious duration.
Agents transition from Recovered back to Susceptible upon waning immunity after the waning duration.
"""

from laser.generic.components import InfectiousIRS
from laser.generic.components import RecoveredRS
from laser.generic.components import Susceptible
from laser.generic.components import TransmissionSI
from laser.generic.shared import State

__all__ = ["InfectiousIRS", "RecoveredRS", "State", "Susceptible", "TransmissionSI"]
