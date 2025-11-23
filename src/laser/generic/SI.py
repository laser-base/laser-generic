"""
Components for an SI model.

Agents transition from Susceptible to Infectious upon infection.
Agents remain in the Infectious state indefinitely (no recovery).
"""

from .components import InfectiousSI as Infectious
from .components import Susceptible
from .components import TransmissionSIX as Transmission
from .shared import State

__all__ = ["Infectious", "State", "Susceptible", "Transmission"]
