"""
Regression test for issue #20: spatial transmission must apply `model.network`
in the correct direction.

The migration step in the ``Transmission*.step`` methods does::

    transfer = ft[:, None] * self.model.network   # transfer[i, j] = ft[i] * net[i, j]
    ft += transfer.sum(axis=0)                    # FOI received by j
    ft -= transfer.sum(axis=1)                    # FOI sent from i

An earlier version used ``transfer = ft * self.model.network`` which (by numpy
broadcasting) gave ``transfer[i, j] = ft[j] * net[i, j]`` — sending infectivity
in the wrong direction for asymmetric networks.

To distinguish the two, this test installs a strictly one-way network
``net = [[0, w], [0, 0]]`` and seeds infectious agents only in node 0. Under the
correct formulation, infectivity flows 0 → 1; under the buggy one, ``transfer``
is all zeros because node 1 has no infectivity to "send back".

The same migration block lives in ``TransmissionSIx`` (SI), ``TransmissionSI``
(SIS/SIR/SIRS), and ``TransmissionSE`` (SEIR/SEIRS). All three are covered.
"""

import unittest
from pathlib import Path

import laser.core.distributions as dists
import numpy as np
from laser.core import PropertySet
from laser.core.random import seed as set_seed

from laser.generic import Model
from laser.generic import SI
from laser.generic.components import (
    Exposed,
    InfectiousIR,
    Recovered,
    Susceptible,
    TransmissionSE,
    TransmissionSI,
    TransmissionSIx,
)

try:
    from tests.utils import stdgrid
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from utils import stdgrid


POPULATION = 10_000
BETA = 0.3
TRANSFER_WEIGHT = 0.1
INF_MEAN = 7.0
EXP_MEAN = 5.0
SEED = 42


def _two_node_scenario():
    """1×2 grid, node 0 fully infectious, node 1 fully susceptible."""
    scenario = stdgrid(M=1, N=2, population_fn=lambda r, c: POPULATION)
    scenario["S"] = [0, POPULATION]
    scenario["I"] = [POPULATION, 0]
    return scenario


def _make_model(scenario, network):
    params = PropertySet({"nticks": 1, "beta": BETA})
    model = Model(scenario, params)
    # Replace the gravity-derived network with our strictly directional one.
    model.network = network
    return model


class _MigrationDirectionMixin:
    """
    Shared assertions for the three transmission components. Subclasses must
    define ``_attach_components(model)`` to wire the right disease components.
    """

    def _attach_components(self, model):  # pragma: no cover - overridden
        raise NotImplementedError

    def _run_with_network(self, network):
        model = _make_model(_two_node_scenario(), network)
        self._attach_components(model)
        model.run("issue #20 migration direction")
        return model.nodes.forces[0]

    def test_forward_one_way_network(self):
        """
        network[0, 1] = w, network[1, 0] = 0. With all infectivity in node 0,
        forces[1] must receive β·w and forces[0] must retain β·(1 − w).
        Under the buggy code transfer is all zeros and forces[1] stays 0.
        """
        network = np.array([[0.0, TRANSFER_WEIGHT], [0.0, 0.0]], dtype=np.float32)
        forces = self._run_with_network(network)
        bare_foi = BETA  # I[0]/N == 1 in node 0 before migration

        assert forces[1] > 0.0, (
            f"Infectivity must flow from node 0 to node 1, but forces[1]={forces[1]}."
            " This is the failure mode of issue #20."
        )
        np.testing.assert_allclose(forces[1], bare_foi * TRANSFER_WEIGHT, rtol=1e-5)
        np.testing.assert_allclose(forces[0], bare_foi * (1.0 - TRANSFER_WEIGHT), rtol=1e-5)

    def test_reverse_one_way_network(self):
        """
        Mirror case: network[1, 0] = w only. Node 0 has no outbound weight, so its
        FOI must stay at β and node 1 must stay at 0. Under the buggy code,
        ``transfer[1, 0] = ft[0] * w`` would be nonzero and forces[0] would leak.
        """
        network = np.array([[0.0, 0.0], [TRANSFER_WEIGHT, 0.0]], dtype=np.float32)
        forces = self._run_with_network(network)
        bare_foi = BETA

        np.testing.assert_allclose(forces[0], bare_foi, rtol=1e-5)
        np.testing.assert_allclose(forces[1], 0.0, atol=1e-7)


class TestTransmissionSIxMigration(_MigrationDirectionMixin, unittest.TestCase):
    """SI model: ``TransmissionSIx``."""

    def setUp(self):
        set_seed(SEED)

    def _attach_components(self, model):
        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=None),
        ]


class TestTransmissionSIMigration(_MigrationDirectionMixin, unittest.TestCase):
    """SIS/SIR/SIRS model: ``TransmissionSI``."""

    def setUp(self):
        set_seed(SEED)

    def _attach_components(self, model):
        infdurdist = dists.normal(loc=INF_MEAN, scale=2.0)
        model.components = [
            Susceptible(model),
            InfectiousIR(model, infdurdist),
            Recovered(model),
            TransmissionSI(model, infdurdist),
        ]


class TestTransmissionSEMigration(_MigrationDirectionMixin, unittest.TestCase):
    """SEIR/SEIRS model: ``TransmissionSE``."""

    def setUp(self):
        set_seed(SEED)

    def _attach_components(self, model):
        scenario = model.scenario
        scenario["E"] = [0, 0]
        scenario["R"] = [0, 0]
        expdurdist = dists.normal(loc=EXP_MEAN, scale=1.0)
        infdurdist = dists.normal(loc=INF_MEAN, scale=2.0)
        model.components = [
            Susceptible(model),
            Exposed(model, expdurdist, infdurdist),
            InfectiousIR(model, infdurdist),
            Recovered(model),
            TransmissionSE(model, expdurdist),
        ]


if __name__ == "__main__":
    unittest.main()
