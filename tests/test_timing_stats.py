"""
Regression test for issue #22: ``TimingStats`` total elapsed should track wall
clock, not 2× wall clock.

The reporter saw a 10-minute wall-clock run report ~20 minutes of total elapsed
time. The current ``_TimingStats`` design tracks inclusive vs exclusive time
per node — a node's inclusive elapsed should equal the wrapped wall-clock time
for that block. If the hierarchical accounting double-counts somewhere, this
test will catch it.

The check compares the elapsed time at the ``"Running Simulation: ..."`` node
(which wraps the entire ``model.run()`` tick loop) against ``time.perf_counter``
around the same call.
"""

import time
import unittest
import uuid
from pathlib import Path

from laser.core import PropertySet
from laser.core.random import seed as set_seed

from laser.generic import Model, SI
from laser.generic.components import TransmissionSIx
from laser.generic.utils import TimingStats

try:
    from tests.utils import stdgrid
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from utils import stdgrid


POPULATION = 20_000
NTICKS = 100
SEED = 42


def _find_node(root, label):
    """Walk the TimingStats tree looking for a node whose label matches."""
    if root.label == label:
        return root
    for child in root.children.values():
        hit = _find_node(child, label)
        if hit is not None:
            return hit
    return None


class TestTimingStatsTracksWallClock(unittest.TestCase):
    """The total elapsed reported by TimingStats must not double-count."""

    def setUp(self):
        set_seed(SEED)

    def test_run_elapsed_within_tolerance_of_wall_clock(self):
        # Unique label so we can find our node in the singleton TimingStats tree.
        unique = f"#22 timing check {uuid.uuid4().hex[:8]}"

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = POPULATION - 100
        scenario["I"] = 100

        params = PropertySet({"nticks": NTICKS, "beta": 0.3, "prng_seed": SEED})
        model = Model(scenario, params)
        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=None),
        ]

        # Warm up numba kernels so JIT cost doesn't dominate the comparison.
        model.run(unique + " warmup")

        # Real measurement run.
        wall_start = time.perf_counter_ns()
        model.run(unique)
        wall_elapsed = time.perf_counter_ns() - wall_start

        node = _find_node(TimingStats.root, f"Running Simulation: {unique}")
        assert node is not None, "Couldn't locate the run's TimingStats node — test plumbing issue."

        reported = node.elapsed
        ratio = reported / wall_elapsed
        # Wall clock and reported elapsed should be approximately equal. Allow
        # generous tolerance for Python/numba/scheduler noise: 0.7 ≤ ratio ≤ 1.3.
        # Issue #22's symptom was ratio ≈ 2; a sane fix puts it well inside this band.
        assert 0.7 <= ratio <= 1.3, (
            f"TimingStats elapsed for run ({reported / 1e6:.1f} ms) diverges from "
            f"wall clock ({wall_elapsed / 1e6:.1f} ms): ratio={ratio:.2f}. "
            "Issue #22 symptom is ratio ≈ 2."
        )


if __name__ == "__main__":
    unittest.main()
