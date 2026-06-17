"""
Reproducibility test for issue #55: two runs with the same ``seed`` must
produce identical state.

``Model.__init__`` now calls ``set_seed(prng_seed)`` with the seed pulled from
``params.prng_seed`` / ``prngseed`` / ``seed`` (in that order). If the seed
mechanism ever breaks again, this test catches it before it ships.
"""

import unittest
from pathlib import Path

import numpy as np
from laser.core import PropertySet

from laser.generic import Model, SI
from laser.generic.components import TransmissionSIx

try:
    from tests.utils import stdgrid
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from utils import stdgrid


POPULATION = 10_000
INITIAL_INFECTED = 50
NTICKS = 50
BETA = 0.3
SEED = 12345


def _run_si_once(seed_key):
    """Build and run an SI model with the seed under the given param key."""
    scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
    scenario["S"] = POPULATION - INITIAL_INFECTED
    scenario["I"] = INITIAL_INFECTED

    params = PropertySet({"nticks": NTICKS, "beta": BETA, seed_key: SEED})
    model = Model(scenario, params)
    model.components = [
        SI.Susceptible(model),
        SI.Infectious(model),
        TransmissionSIx(model, seasonality=None),
    ]
    model.run(f"reproducibility[{seed_key}]")
    return model.nodes.S.copy(), model.nodes.I.copy(), model.people.state.copy()


class TestSameSeedSameResult(unittest.TestCase):
    """Two runs with the same ``seed`` parameter must produce identical state."""

    def test_seed_param_makes_runs_reproducible(self):
        S1, I1, state1 = _run_si_once("seed")
        S2, I2, state2 = _run_si_once("seed")

        np.testing.assert_array_equal(S1, S2, err_msg="nodes.S diverged between runs with same seed (issue #55).")
        np.testing.assert_array_equal(I1, I2, err_msg="nodes.I diverged between runs with same seed (issue #55).")
        np.testing.assert_array_equal(state1, state2, err_msg="people.state diverged between runs with same seed (issue #55).")

    def test_prng_seed_param_makes_runs_reproducible(self):
        # Model.__init__ also accepts prng_seed; cover that alias explicitly.
        S1, I1, _ = _run_si_once("prng_seed")
        S2, I2, _ = _run_si_once("prng_seed")
        np.testing.assert_array_equal(S1, S2)
        np.testing.assert_array_equal(I1, I2)


if __name__ == "__main__":
    unittest.main()
