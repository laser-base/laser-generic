__version__ = "1.1.0"

from laser.generic.immunization import ImmunizationCampaignEx as ImmunizationCampaign
from laser.generic.immunization import RoutineImmunizationEx as RoutineImmunization
from laser.generic.importation import Infect_Random_Agents
from laser.generic.model import Model
from laser.generic.shared import State

from . import SEIR
from . import SEIRS
from . import SI
from . import SIR
from . import SIRS
from . import SIS

__all__ = [
    "SEIR",
    "SEIRS",
    "SI",
    "SIR",
    "SIRS",
    "SIS",
    "ImmunizationCampaign",
    "Infect_Random_Agents",
    "Model",
    "RoutineImmunization",
    "State",
]
