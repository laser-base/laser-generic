"""
Export required components for an SIRS model.

Agents transition from Susceptible to Infectious upon infection.
Agents transition from Infectious to Recovered upon recovery after the infectious duration.
Agents transition from Recovered back to Susceptible upon waning immunity after the waning duration.
"""

from laser.generic.components import InfectiousIRS as Infectious
from laser.generic.components import RecoveredRS as Recovered
from laser.generic.components import Susceptible
from laser.generic.components import TransmissionSI as Transmission
from laser.generic.shared import State

__all__ = ["Infectious", "Recovered", "State", "Susceptible", "Transmission"]
