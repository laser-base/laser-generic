"""
Export required components for an SEIR model.

Agents transition from Susceptible to Exposed upon infection, with an incubation duration.
Agents transition from Exposed to Infectious after the incubation period and are infectious for a duration.
Agents transition from Infectious to Recovered after the infectious period.
Agents remain in the Recovered state indefinitely (no waning immunity).
"""

from laser.generic.components import Exposed
from laser.generic.components import InfectiousIR as Infectious
from laser.generic.components import Recovered
from laser.generic.components import Susceptible
from laser.generic.components import TransmissionSE as Transmission
from laser.generic.shared import State

__all__ = ["Exposed", "Infectious", "Recovered", "State", "Susceptible", "Transmission"]
