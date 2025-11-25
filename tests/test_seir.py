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

from laser.generic import SEIR
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
EXPOSED_DURATION_SHAPE = 4.5
EXPOSED_DURATION_SCALE = 1.0
INFECTIOUS_DURATION_MEAN = 7.0


def build_model(m, n, pop_fn, init_infected=0, init_recovered=0, birthrates=None, pyramid=None, survival=None):
    """
    Helper function: build a complete SEIR model with configurable demography.
    Creates Susceptible, Exposed, Infectious, and Recovered components plus
    optional birth and mortality processes.
    """
    scenario = stdgrid(M=m, N=n, population_fn=pop_fn)
    scenario["S"] = scenario["population"]
    scenario["E"] = 0
    scenario["S"] -= init_infected
    scenario["I"] = init_infected
    scenario["S"] -= init_recovered
    scenario["R"] = init_recovered

    beta = R0 / INFECTIOUS_DURATION_MEAN
    params = PropertySet({"nticks": NTICKS, "beta": beta})

    with ts.start("Model Initialization"):
        model = Model(scenario, params, birthrates=birthrates)

        expdist = dists.gamma(shape=EXPOSED_DURATION_SHAPE, scale=EXPOSED_DURATION_SCALE)
        infdist = dists.normal(loc=INFECTIOUS_DURATION_MEAN, scale=2)

        s = SEIR.Susceptible(model)
        e = SEIR.Exposed(model, expdist, infdist)
        i = SEIR.Infectious(model, infdist)
        r = SEIR.Recovered(model)
        tx = SEIR.Transmission(model, expdist)

        if birthrates is not None:
            assert pyramid is not None, "Pyramid must be provided for vital dynamics."
            assert survival is not None, "Survival function must be provided for vital dynamics."
            births = BirthsByCBR(model, birthrates, pyramid)
            mortality = MortalityByEstimator(model, survival)
            model.components = [s, e, i, r, tx, births, mortality]
        else:
            model.components = [s, e, i, r, tx]

        model.validating = VALIDATING

    return model


class Default(unittest.TestCase):
    def test_single(self):
        """
        Feature: Single-node deterministic SEIR model
        --------------------------------------------------
        Validates:
          • Infection latency (E state) prior to infectiousness.
          • Proper sequencing of transitions S→E→I→R.
          • Final recovered fraction consistent with R₀ = 1.386.
          • Population mass conservation (S+E+I+R constant).

        Configuration:
          Nodes: 1
          Population: 100,000
          Initial infections: 10
          Exposure duration: Gamma(shape=4.5, scale=1)
          Infectious duration: Normal(mean=7, sd=2)
          Simulation length: 365 ticks

        Expected Outcomes / Invariants:
          • Infection curve rises and decays with clear latent period.
          • E and I trajectories overlap correctly (delayed onset).
          • Population constant to within 0.01%.
          • Final R fraction ≈ 0.5 ± 0.05.
        """
        with ts.start("test_single_node"):
            model = build_model(1, 1, lambda x, y: 100_000, init_infected=10)
            model.run("SEIR Single Node")

            I_series = model.nodes.I.sum(axis=1)
            E_series = model.nodes.E.sum(axis=1)
            # R_series = model.nodes.R.sum(axis=1)

            # Quantitative checks
            assert E_series.max() > 0, "No exposed cases observed."
            assert np.argmax(E_series) < np.argmax(I_series), "E should peak before I."
            peak_I = np.argmax(I_series)
            assert I_series[-1] < I_series[peak_I] * 0.5, "I should decline post-peak."

            N0 = (model.nodes.S[0] + model.nodes.E[0] + model.nodes.I[0] + model.nodes.R[0]).sum()
            NT = (model.nodes.S[-1] + model.nodes.E[-1] + model.nodes.I[-1] + model.nodes.R[-1]).sum()
            assert abs(NT - N0) / N0 < 1e-4, "Population not conserved (ΔN>0.01%)."

            final_R_frac = model.nodes.R[-1].sum() / N0
            assert 0.45 <= final_R_frac <= 0.55, f"Final attack fraction {final_R_frac:.3f} out of expected range."

    def test_grid(self):
        """
        Feature: Spatial 2-D SEIR model with births and deaths
        --------------------------------------------------
        Validates:
          • Spatial epidemic propagation with exposed delay dynamics.
          • Integration of birth and mortality processes.
          • Stability of total population under demographic turnover.
          • Infection prevalence within epidemiologically realistic bounds.

        Configuration:
          Grid: 10x10 nodes (100 total)
          Population: 10,000-1,000,000 per node
          Exposure: Gamma(shape=4.5, scale=1)
          Infectious: Normal(mean=7, sd=2)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • Mean prevalence (I/N) ≤ 0.5.
          • Population drift ≤ ±10%.
          • E precedes I temporally.
          • No negative state counts.
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
            model.run("SEIR Grid")

            I_series = model.nodes.I.sum(axis=1)
            E_series = model.nodes.E.sum(axis=1)
            N_series = (model.nodes.S + model.nodes.E + model.nodes.I + model.nodes.R).sum(axis=1)
            mean_prev = (model.nodes.I / (model.nodes.S + model.nodes.E + model.nodes.I + model.nodes.R + 1e-9)).mean()
            pop_change = (N_series[-1] - N_series[0]) / N_series[0]

            assert np.all(model.nodes.S >= 0)
            assert np.all(model.nodes.E >= 0)
            assert np.all(model.nodes.I >= 0)
            assert np.all(model.nodes.R >= 0)
            assert mean_prev <= 0.5, f"Mean prevalence {mean_prev:.3f} > 0.5"
            assert abs(pop_change) < 0.1, f"Population drift {pop_change * 100:.2f}% >10%"
            assert np.argmax(E_series) < np.argmax(I_series), "E should peak before I."

    def test_linear(self):
        """
        Feature: One-dimensional (linear) SEIR model
        --------------------------------------------------
        Validates:
          • Epidemic spread along a 1xN chain.
          • Latent exposure delays and sequential infection wave.
          • Population stability and bounded infection.

        Configuration:
          Layout: 1x10 chain
          Exposure: Gamma(shape=4.5, scale=1)
          Infectious: Normal(mean=7, sd=2)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • Infection grows and then decays (non-static dynamics).
          • Total population drift <5%.
          • E precedes I across chain nodes.
          • All counts non-negative.
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
            model.run("SEIR Linear")

            I_series = model.nodes.I.sum(axis=1)
            E_series = model.nodes.E.sum(axis=1)
            pop_series = (model.nodes.S + model.nodes.E + model.nodes.I + model.nodes.R).sum(axis=1)
            pop_change = (pop_series[-1] - pop_series[0]) / pop_series[0]

            assert np.all(model.nodes.S >= 0)
            assert np.all(model.nodes.E >= 0)
            assert np.all(model.nodes.I >= 0)
            assert np.all(model.nodes.R >= 0)
            assert I_series.max() > I_series[0] * 1.5, "Epidemic growth too weak."
            assert I_series[-1] < I_series.max() * 0.8, "Epidemic did not decline."
            assert abs(pop_change) < 0.05, f"Population drift {pop_change * 100:.2f}% >5%"
            assert np.argmax(E_series) < np.argmax(I_series), "E should peak before I."


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
