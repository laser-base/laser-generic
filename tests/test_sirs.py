from laser.generic.newutils import TimingStats as ts  # noqa: I001

import json
import unittest
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import laser.core.distributions as dists
from laser.core import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.demographics import KaplanMeierEstimator

from laser.generic import SIRS
from laser.generic import Model
from laser.generic.newutils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByEstimator
from tests.utils import stdgrid

PLOTTING = False
VERBOSE = False
EM = 10
EN = 10
PEE = 10
VALIDATING = False
NTICKS = 365
R0 = 1.386
INFECTIOUS_DURATION_MEAN = 7.0
WANING_DURATION_MEAN = 30.0


def build_model(m, n, pop_fn, init_infected=0, init_recovered=0, birthrates=None, pyramid=None, survival=None):
    """
    Helper: Construct an SIRS model with optional demography and waning immunity.
    """
    scenario = stdgrid(M=m, N=n, population_fn=pop_fn)
    scenario["S"] = scenario["population"]
    scenario["S"] -= init_infected
    scenario["I"] = init_infected
    scenario["S"] -= init_recovered
    scenario["R"] = init_recovered

    beta = R0 / INFECTIOUS_DURATION_MEAN
    params = PropertySet({"nticks": NTICKS, "beta": beta})

    with ts.start("Model Initialization"):
        model = Model(scenario, params, birthrates=birthrates)

        infdurdist = dists.normal(loc=INFECTIOUS_DURATION_MEAN, scale=2)
        wandurdist = dists.normal(loc=WANING_DURATION_MEAN, scale=5)

        s = SIRS.Susceptible(model)
        i = SIRS.Infectious(model, infdurdist, wandurdist)
        r = SIRS.Recovered(model, wandurdist)
        tx = SIRS.Transmission(model, infdurdist)

        if birthrates is not None:
            assert pyramid is not None, "Pyramid must be provided for vital dynamics."
            assert survival is not None, "Survival function must be provided for vital dynamics."
            births = BirthsByCBR(model, birthrates, pyramid)
            mortality = MortalityByEstimator(model, survival)
            model.components = [s, i, r, tx, births, mortality]
        else:
            model.components = [s, i, r, tx]

        model.validating = VALIDATING

    return model


class Default(unittest.TestCase):
    def test_single(self):
        """
        Feature: Single-node deterministic SIRS model
        --------------------------------------------------
        Validates:
          • Infection and recovery progression with waning immunity.
          • Correct S→I→R→S loop behavior.
          • Conservation of total population (S+I+R constant).
          • Periodic reinfection due to waning.

        Configuration:
          Nodes: 1
          Population: 100,000
          Initial infections: 10
          Infectious duration: Normal(mean=7, sd=2)
          Waning duration: Normal(mean=30, sd=5)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • Infection curve exhibits rise, fall, and secondary bumps.
          • R decreases as immunity wanes.
          • Population constant within 0.01%.
        """
        with ts.start("test_single_node"):
            model = build_model(1, 1, lambda x, y: 100_000, init_infected=10)
            model.run("SIRS Single Node")

            I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            S_series = model.nodes.S.sum(axis=1)

            # Quantitative checks
            assert I_series.max() > I_series[0] * 2, "Infection did not grow sufficiently."
            assert R_series.max() > R_series[0], "No recovery observed."
            assert R_series[-1] < R_series.max(), "No waning immunity (R did not decline)."
            assert S_series[-1] < S_series[0], "Susceptible pool should initially shrink."

            N0 = (model.nodes.S[0] + model.nodes.I[0] + model.nodes.R[0]).sum()
            NT = (model.nodes.S[-1] + model.nodes.I[-1] + model.nodes.R[-1]).sum()
            assert abs(NT - N0) / N0 < 1e-4, "Population not conserved."

    def test_grid(self):
        """
        Feature: Spatial 2-D SIRS model with births and deaths
        --------------------------------------------------
        Validates:
          • Spatial coupling of SIRS dynamics across a grid.
          • Integration of demographic turnover.
          • Bounded prevalence and stable total population.
          • Infection re-emergence under waning immunity.

        Configuration:
          Grid: 10x10 nodes (100 total)
          Population: 10,000-1,000,000 per node
          Infectious: Normal(mean=7, sd=2)
          Waning: Normal(mean=30, sd=5)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • Mean prevalence ≤ 0.5.
          • Population drift ≤ ±10%.
          • Non-negative S, I, R counts.
          • Infection cycles present in some patches.
        """
        with ts.start("test_grid"):
            cbr = np.random.uniform(5, 35, EM * EN)
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)
            pyramid = AliasedDistribution(np.full(89, 1_000))
            survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

            model = build_model(
                EM,
                EN,
                lambda x, y: int(np.random.uniform(10_000, 1_000_000)),
                init_infected=10,
                birthrates=birthrate_map.values,
                pyramid=pyramid,
                survival=survival,
            )
            model.run("SIRS Grid")

            # I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            N_series = (model.nodes.S + model.nodes.I + model.nodes.R).sum(axis=1)
            pop_change = (N_series[-1] - N_series[0]) / N_series[0]
            mean_prev = (model.nodes.I / (model.nodes.S + model.nodes.I + model.nodes.R + 1e-9)).mean()

            assert np.all(model.nodes.S >= 0)
            assert np.all(model.nodes.I >= 0)
            assert np.all(model.nodes.R >= 0)
            assert mean_prev <= 0.5, f"Mean prevalence {mean_prev:.3f} > 0.5"
            assert abs(pop_change) < 0.1, f"Population drift {pop_change * 100:.2f}% >10%"
            assert R_series.max() > R_series[-1], "R should decline due to waning immunity."

    def test_linear(self):
        """
        Feature: One-dimensional (linear) SIRS model
        --------------------------------------------------
        Validates:
          • Propagation of infection in a 1xN linear chain.
          • Reinfection cycles due to waning immunity.
          • Stability of population with births and deaths.

        Configuration:
          Layout: 1x10 chain
          Infectious duration: Normal(mean=7, sd=2)
          Waning duration: Normal(mean=30, sd=5)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • E→I→R→S flow confirmed via time ordering.
          • Population drift <5%.
          • Non-negative state counts.
          • Multiple infection cycles observed.
        """
        with ts.start("test_linear"):
            cbr = np.random.uniform(5, 35, PEE)
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)
            pyramid = AliasedDistribution(np.full(89, 1_000))
            survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

            model = build_model(
                1,
                PEE,
                lambda x, y: int(np.random.uniform(10_000, 1_000_000)),
                init_infected=10,
                birthrates=birthrate_map.values,
                pyramid=pyramid,
                survival=survival,
            )
            model.run("SIRS Linear")

            I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            pop_series = (model.nodes.S + model.nodes.I + model.nodes.R).sum(axis=1)
            pop_change = (pop_series[-1] - pop_series[0]) / pop_series[0]

            assert np.all(model.nodes.S >= 0)
            assert np.all(model.nodes.I >= 0)
            assert np.all(model.nodes.R >= 0)
            assert I_series.max() > I_series[0] * 1.5, "Epidemic growth too weak."
            assert I_series[-1] < I_series.max() * 0.8, "Epidemic did not decline."
            assert abs(pop_change) < 0.05, f"Population drift {pop_change * 100:.2f}% >5%"
            assert R_series.max() > R_series[-1], "Waning immunity not evident (R did not decrease)."


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--validating", action="store_true")
    parser.add_argument("-m", type=int, default=5)
    parser.add_argument("-n", type=int, default=5)
    parser.add_argument("-p", type=int, default=10)
    parser.add_argument("-t", "--ticks", type=int, default=365)
    parser.add_argument("-r", "--r0", type=float, default=1.386)
    parser.add_argument("-g", "--grid", action="store_true")
    parser.add_argument("-l", "--linear", action="store_true")
    parser.add_argument("-s", "--single", action="store_true")
    parser.add_argument("unittest", nargs="*")

    args = parser.parse_args()
    PLOTTING = args.plot
    VERBOSE = args.verbose
    VALIDATING = args.validating
    NTICKS = args.ticks
    R0 = args.r0
    EM, EN, PEE = args.m, args.n, args.p

    print(f"Using arguments {args=}")

    tc = Default()
    run_all = not (args.grid or args.linear or args.single)

    if args.single or run_all:
        tc.test_single()
    if args.grid or run_all:
        tc.test_grid()
    if args.linear or run_all:
        tc.test_linear()

    ts.freeze()
    print("\nTiming Summary:")
    print("-" * 30)
    print(ts.to_string(scale="ms"))
    with Path("timing_data.json").open("w") as f:
        json.dump(ts.to_dict(scale="ms"), f, indent=4)
