from enum import Enum

import numpy as np
from laser.core.demographics import AliasedDistribution
from laser.core.demographics import KaplanMeierEstimator

__all__ = ["State", "sample_dobs", "sample_dods"]


def sample_dobs(pyramid: AliasedDistribution, dobs: np.ndarray, tick: int = 0) -> None:
    # Get years of age sampled from the population pyramid
    dobs[:] = pyramid.sample(len(dobs)).astype(dobs.dtype)
    dobs *= 365  # Convert years to days
    dobs += np.random.randint(0, 365, size=len(dobs))  # add some noise within the year
    # pyramid.sample actually returned ages. Turn them into dobs by treating them
    # as days before today.
    dobs[:] = tick - dobs

    return


def sample_dods(dobs: np.ndarray, survival: KaplanMeierEstimator, tick: int, dods: np.ndarray) -> None:
    # An agent's age is (tick - dob).
    dods[:] = survival.predict_age_at_death(tick - dobs).astype(dods.dtype)
    dods += dobs  # Convert days until death to a date of death.

    return


class State(Enum):
    DECEASED = -1
    SUSCEPTIBLE = 0
    EXPOSED = 1
    INFECTIOUS = 2
    RECOVERED = 3

    def __new__(cls, value):
        obj = object.__new__(cls)
        obj._value_ = np.int8(value)
        return obj
