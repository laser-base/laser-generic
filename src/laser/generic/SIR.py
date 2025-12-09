"""
Export required components for an SIR model.

Agents transition from Susceptible to Infectious upon infection and are infectious for a duration.
Agents transition from Infectious to Recovered upon recovery.
Agents remain in the Recovered state indefinitely (no waning immunity).
"""

from laser.generic.components import InfectiousIR
from laser.generic.components import Recovered
from laser.generic.components import Susceptible
from laser.generic.components import TransmissionSI
from laser.generic.shared import State

__all__ = ["InfectiousIR", "Recovered", "State", "Susceptible", "TransmissionSI"]
