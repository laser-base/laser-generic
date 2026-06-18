"""
Capacity-on-late-births test for issue #114.

When ``Model`` is built without a ``birthrates`` argument, capacity is computed
as if the population were static — there is no headroom for new agents. If a
``BirthsByCBR`` component is then attached and run, the LaserFrame runs out of
space the first time births fire.

#114 proposes that ``Model`` either auto-detect vital-dynamics components and
adjust capacity, or fail with a clear error before ``run()`` starts. Either
outcome would make this test pass; today it xfails because the LaserFrame
overflow happens mid-tick during the run.
"""

import unittest
from pathlib import Path

import numpy as np
import pytest
from laser.core import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.random import seed as set_seed

from laser.generic import Model, SI
from laser.generic.components import TransmissionSIx
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR

try:
    from tests.utils import stdgrid
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from utils import stdgrid


POPULATION = 5_000
NTICKS = 365
CBR = 30.0  # crude birth rate per 1000 per year — high so births fire quickly
SEED = 42


class TestModelCapacityWithLateAddedBirths(unittest.TestCase):
    """Model built without birthrates must accommodate a later-added BirthsByCBR."""

    def setUp(self):
        set_seed(SEED)

    def test_capacity_matches_population_when_no_birthrates_given(self):
        """Smoking-gun assertion: with no birthrates, capacity equals initial population."""
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = POPULATION - 10
        scenario["I"] = 10

        params = PropertySet({"nticks": NTICKS, "beta": 0.0, "prng_seed": SEED})
        model = Model(scenario, params)

        # This is the precondition that makes #114 a real problem. If Model
        # learns to allocate headroom proactively, this assertion changes.
        assert model.people.capacity == POPULATION, (
            f"With no birthrates, Model allocates exactly the initial population "
            f"({POPULATION}) — leaving no room for later-added BirthsByCBR. "
            "This is the precondition described in issue #114."
        )

    @pytest.mark.xfail(
        reason=(
            "Issue #114: Model built without birthrates allocates capacity = "
            "initial population, so a BirthsByCBR component attached afterward "
            "overflows the LaserFrame the first tick births fire — currently "
            "raises ValueError('frame.add() exceeds capacity ...') mid-run. "
            "The fix is either to auto-detect vital-dynamics components and "
            "grow capacity, or to raise a clear error at run() start. Remove "
            "this xfail when either is implemented."
        ),
        raises=ValueError,
        strict=True,
    )
    def test_model_runs_when_births_attached_after_construction(self):
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = POPULATION - 10
        scenario["I"] = 10

        params = PropertySet({"nticks": NTICKS, "beta": 0.0, "prng_seed": SEED})

        # Construct Model with no birthrates — capacity will be tight.
        model = Model(scenario, params)

        # Build a positive birthrate map and a trivial age pyramid AFTER the
        # Model exists, mimicking the user error pattern #114 describes.
        cbr_map = ValuesMap.from_scalar(CBR, NTICKS, model.nodes.count).values
        pyramid = AliasedDistribution(np.full(89, 1_000))

        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=None),
            BirthsByCBR(model, birthrates=cbr_map, pyramid=pyramid),
        ]

        model.run("Issue #114: late-added births")


if __name__ == "__main__":
    unittest.main()
