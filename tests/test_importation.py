"""
Construction and accounting tests for the importation components.

These pair with the bug filing that ``importation.py`` references the stale
``model.patches.*`` namespace (the LaserFrame is now ``model.nodes``). Both
tests xfail on ``main``; they should flip to passing when the namespace is
swept.

What the bug does today:

- ``Infect_Agents_In_Patch.__init__`` falls back to ``np.arange(model.patches.count)``
  when no ``importation_patchlist`` is supplied — raises ``AttributeError``
  immediately.
- ``Infect_Random_Agents.__call__`` guards on ``hasattr(model.patches, "cases_test")``
  which silently returns False, so the ``cases_test`` / ``susceptibility_test``
  accounting branch never runs even when those channels are installed.
"""

import unittest
from pathlib import Path

import numpy as np
import pytest
from laser.core import PropertySet
from laser.core.random import seed as set_seed

from laser.generic import Model, SI
from laser.generic.components import TransmissionSIx
from laser.generic.importation import Infect_Agents_In_Patch, Infect_Random_Agents

try:
    from tests.utils import stdgrid
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from utils import stdgrid


POPULATION = 5_000
NTICKS = 20
SEED = 42


def _build_si_model(num_nodes=2, nticks=NTICKS):
    scenario = stdgrid(M=1, N=num_nodes, population_fn=lambda r, c: POPULATION)
    scenario["S"] = POPULATION
    scenario["I"] = 0
    params = PropertySet(
        {
            "nticks": nticks,
            "beta": 0.0,  # disable native transmission so importation is the only source
            "importation_period": 5,
            "importation_count": 3,
            "prng_seed": SEED,
        }
    )
    model = Model(scenario, params)
    model.components = [
        SI.Susceptible(model),
        SI.Infectious(model),
        TransmissionSIx(model, seasonality=None),
    ]
    return model


class TestInfectAgentsInPatchConstruction(unittest.TestCase):
    """``Infect_Agents_In_Patch`` must initialize on the current ``Model``."""

    def setUp(self):
        set_seed(SEED)

    @pytest.mark.xfail(
        reason=(
            "importation.py:149 references model.patches.count when "
            "importation_patchlist is not supplied. model.patches no longer "
            "exists on Model — should be model.nodes.count. Filed as the "
            "'stale model.patches.* namespace' bug."
        ),
        raises=AttributeError,
        strict=True,
    )
    def test_init_without_importation_patchlist_defaults_to_all_nodes(self):
        model = _build_si_model(num_nodes=3)
        component = Infect_Agents_In_Patch(model)
        np.testing.assert_array_equal(component.patchlist, np.arange(model.nodes.count))


class TestInfectRandomAgentsCasesTestBranch(unittest.TestCase):
    """The ``cases_test`` / ``susceptibility_test`` accounting branch must trigger
    when those channels are present."""

    def setUp(self):
        set_seed(SEED)

    @pytest.mark.xfail(
        reason=(
            "importation.py:82 guards on hasattr(model.patches, 'cases_test'). "
            "model.patches doesn't exist, so the branch is silently dead even "
            "when cases_test/susceptibility_test are installed on model.nodes. "
            "Filed as the 'stale model.patches.* namespace' bug."
        ),
        strict=True,
    )
    def test_cases_test_channel_updated_on_importation_tick(self):
        model = _build_si_model(num_nodes=2, nticks=NTICKS)
        # Install the census-style channels the importation component expects to update.
        model.nodes.add_vector_property("cases_test", NTICKS + 1, dtype=np.int32)
        model.nodes.add_vector_property("susceptibility_test", NTICKS + 1, dtype=np.int32)
        before_cases = model.nodes.cases_test.copy()

        component = Infect_Random_Agents(model)
        model.components = [*model.components, component]
        model.run("Infect_Random_Agents accounting")

        # Importation fires at ticks 0, 5, 10, 15 with count=3 → 12 total infections.
        delta = model.nodes.cases_test - before_cases
        assert delta.sum() > 0, (
            "cases_test accounting branch never ran; importation.py is gated on "
            "the stale model.patches.* namespace."
        )


if __name__ == "__main__":
    unittest.main()
