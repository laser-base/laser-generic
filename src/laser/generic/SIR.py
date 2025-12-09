"""
Export required components for an SIR model.

Agents transition from Susceptible to Infectious upon infection and are infectious for a duration.
Agents transition from Infectious to Recovered upon recovery.
Agents remain in the Recovered state indefinitely (no waning immunity).
"""

from .components import InfectiousIR
from .components import Recovered
from .components import Susceptible
from .components import TransmissionSI
from .shared import State

__all__ = ["InfectiousIR", "Recovered", "State", "Susceptible", "TransmissionSI"]
