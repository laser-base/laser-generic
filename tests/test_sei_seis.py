"""
Smoke tests for the SEI and SEIS configurations (verification for issue #32).

There is no ``laser.generic.SEI`` or ``SEIS`` package — they're hand-wired from
components (see ``docs/tutorials/notebooks/SEI_and_SEIS_implementations.ipynb``).
The bug in #32 reported that adding ``Exposed`` + ``Infectious*`` produced no
transmission. These tests pin down a working recipe so that regressions in
either the wiring or the underlying components show up in CI.
"""

import unittest
from pathlib import Path

import laser.core.distributions as dists
import numpy as np
import pytest
from laser.core import PropertySet
from laser.core.random import seed as set_seed

from laser.generic import Model
from laser.generic.components import (
    Exposed,
    InfectiousIS,
    InfectiousSI,
    Susceptible,
    TransmissionSE,
)

try:
    from tests.utils import stdgrid
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from utils import stdgrid


POPULATION = 50_000
INITIAL_INFECTED = 500
NTICKS = 200
BETA = 0.4
EXP_MEAN = 5.0
INF_MEAN = 7.0
SEED = 42


def _scenario_with_e_column(states):
    """1x1 grid scenario, pre-populated with the given state columns."""
    scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
    for col, value in states.items():
        scenario[col] = value
    return scenario


class TestSEIModelTransmits(unittest.TestCase):
    """SEI = S → E → I (no recovery). Transmission must actually flow S → E → I."""

    def setUp(self):
        set_seed(SEED)

    @pytest.mark.xfail(
        reason=(
            "Issue #32: Exposed assigns model.people.itimer on E→I transition, but "
            "InfectiousSI never adds the itimer field — the components don't compose. "
            "Either Exposed needs to skip itimer assignment when no recovery component "
            "is present, or an 'InfectiousI' variant that registers itimer without "
            "recovery logic is needed. Remove this xfail when #32 is fixed."
        ),
        raises=AttributeError,
        strict=True,
    )
    def test_sei_transmission_progresses_through_all_states(self):
        scenario = _scenario_with_e_column(
            {
                "S": POPULATION - INITIAL_INFECTED,
                "E": 0,
                "I": INITIAL_INFECTED,
            }
        )

        params = PropertySet({"nticks": NTICKS, "beta": BETA, "prng_seed": SEED})
        model = Model(scenario, params)

        expdurdist = dists.normal(loc=EXP_MEAN, scale=1.0)
        infdurdist = dists.normal(loc=INF_MEAN, scale=2.0)

        # SEI: Exposed handles E→I (and needs infdurdist for assigning itimer).
        # InfectiousSI never recovers (no I→? transition).
        model.components = [
            Susceptible(model),
            Exposed(model, expdurdist, infdurdist),
            InfectiousSI(model),
            TransmissionSE(model, expdurdist),
        ]

        model.run("SEI smoke test")

        S = model.nodes.S
        E = model.nodes.E
        I = model.nodes.I  # noqa: E741

        assert E[NTICKS].sum() > 0 or E[:NTICKS].sum() > 0, (
            "Exposed compartment must have been populated at some point — "
            "transmission appears blocked (issue #32 symptom)."
        )
        assert I[NTICKS].sum() > I[0].sum(), (
            f"Infectious count must grow under SEI: I[0]={I[0].sum()} I[end]={I[NTICKS].sum()}"
        )
        assert S[NTICKS].sum() < S[0].sum(), (
            "Susceptible pool must shrink as agents move S → E → I."
        )

        # Population conservation (no births, no deaths, no recovery exit)
        total = S[NTICKS].sum() + E[NTICKS].sum() + I[NTICKS].sum()
        assert total == POPULATION, (
            f"SEI conserves population: expected {POPULATION}, got {total}"
        )


class TestSEISModelTransmits(unittest.TestCase):
    """SEIS = S → E → I → S (no recovered compartment, agents cycle back to S)."""

    def setUp(self):
        set_seed(SEED)

    def test_seis_transmission_and_return_to_susceptible(self):
        scenario = _scenario_with_e_column(
            {
                "S": POPULATION - INITIAL_INFECTED,
                "E": 0,
                "I": INITIAL_INFECTED,
            }
        )

        params = PropertySet({"nticks": NTICKS, "beta": BETA, "prng_seed": SEED})
        model = Model(scenario, params)

        expdurdist = dists.normal(loc=EXP_MEAN, scale=1.0)
        infdurdist = dists.normal(loc=INF_MEAN, scale=2.0)

        # SEIS: InfectiousIS sends I→S on itimer expiry.
        model.components = [
            Susceptible(model),
            Exposed(model, expdurdist, infdurdist),
            InfectiousIS(model, infdurdist),
            TransmissionSE(model, expdurdist),
        ]

        model.run("SEIS smoke test")

        S = model.nodes.S
        E = model.nodes.E
        I = model.nodes.I  # noqa: E741

        assert E[:NTICKS].sum() > 0, (
            "Exposed compartment must have been populated at some point — "
            "transmission appears blocked (issue #32 symptom)."
        )
        assert I[:NTICKS].sum() > I[0].sum(), (
            "Infectious counts must grow at some point under SEIS — the I→S "
            "loop should sustain transmission longer than under SEI."
        )

        # Population conservation
        total = S[NTICKS].sum() + E[NTICKS].sum() + I[NTICKS].sum()
        assert total == POPULATION, (
            f"SEIS conserves population: expected {POPULATION}, got {total}"
        )


if __name__ == "__main__":
    unittest.main()
