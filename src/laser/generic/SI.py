"""
Components for an SI model.

Agents transition from Susceptible to Infectious upon infection.
Agents remain in the Infectious state indefinitely (no recovery).
"""

from .components import InfectiousSI
from .components import Susceptible
from .components import TransmissionSIX
from .shared import State

__all__ = ["InfectiousSI", "State", "Susceptible", "TransmissionSIX"]
