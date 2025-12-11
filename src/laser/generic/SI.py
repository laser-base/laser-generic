"""
Components for an SI model.

Agents transition from Susceptible to Infectious upon infection.
Agents remain in the Infectious state indefinitely (no recovery).
"""

from laser.generic.components import InfectiousSI as Infectious
from laser.generic.components import Susceptible
from laser.generic.components import TransmissionSIx as Transmission
from laser.generic.shared import State

__all__ = ["Infectious", "State", "Susceptible", "Transmission"]
