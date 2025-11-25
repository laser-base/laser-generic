import json
import sys
import unittest
from argparse import ArgumentParser
from pathlib import Path

import laser.core.distributions as dists
import numpy as np
import pytest
from laser.core import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.demographics import KaplanMeierEstimator

from laser.generic import SIS
from laser.generic import Model
from laser.generic.newutils import TimingStats as ts
from laser.generic.newutils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByEstimator
from utils import base_maps
from utils import stdgrid

PLOTTING = False
VERBOSE = False
EM = 10
EN = 10
PEE = 10
VALIDATING = False
NTICKS = 365


@pytest.mark.feature("spatial_grid_sis_model")
@pytest.mark.feature("demography_births_deaths")
@pytest.mark.feature("infection_transmission_progression")
class Default(unittest.TestCase):
    def test_grid(self):
        """Validate SIS grid dynamics on a 2-D spatial lattice."""
        with ts.start("test_grid"):
            grd = stdgrid(
                M=EM,
                N=EN,
                node_size_km=10,
                population_fn=lambda x, y: int(np.random.uniform(10_000, 1_000_000)),
                origin_x=-119.204167,
                origin_y=40.786944,
            )
            scenario = grd
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10

            # --- Basic population sanity ---
            assert np.all(scenario["S"] >= 0)
            assert np.all(scenario["I"] >= 0)
            np.testing.assert_array_equal(scenario["S"] + scenario["I"], scenario["population"])

            # Birthrates and parameters
            cbr = np.random.uniform(5, 35, len(scenario))  # per 1,000 per year
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)

            R0 = 1.2
            infectious_duration_mean = 7.0
            beta = R0 / infectious_duration_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrates=birthrate_map)

                infdist = dists.normal(loc=infectious_duration_mean, scale=2)
                pyramid = AliasedDistribution(np.full(89, 1_000))
                survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

                s = SIS.Susceptible(model)
                i = SIS.Infectious(model, infdist)
                tx = SIS.Transmission(model, infdist)
                births = BirthsByCBR(model, birthrate_map, pyramid)
                mortality = MortalityByEstimator(model, survival)
                model.components = [s, i, tx, births, mortality]

                model.validating = VALIDATING

            # Run model
            model.run(f"SIS Grid ({model.people.count:,}/{model.nodes.count:,})")

        # --- Post-simulation checks ---
        total_pop = model.nodes.S[-1].sum() + model.nodes.I[-1].sum()
        assert total_pop > 0, "Population should remain positive"
        assert np.all(model.nodes.S >= 0)
        assert np.all(model.nodes.I >= 0)
        # infection prevalence should not exceed 1.0
        assert np.all(model.nodes.I / (model.nodes.I + model.nodes.S + 1e-9) <= 1.0)

        if VERBOSE:
            print(model.people.describe("People"))
            print(model.nodes.describe("Nodes"))

        if PLOTTING:
            ibm = np.random.choice(len(base_maps))
            model.basemap_provider = base_maps[ibm]
            print(f"Using basemap: {model.basemap_provider.name}")
            model.plot()

    @pytest.mark.feature("linear_chain_sis_model")
    def test_linear(self):
        """Validate SIS model on a 1-D linear arrangement of patches."""
        with ts.start("test_linear"):
            lin = stdgrid(
                M=1,
                N=PEE,
                node_size_km=10,
                population_fn=lambda x, y: int(np.random.uniform(10_000, 1_000_000)),
                origin_x=-119.204167,
                origin_y=40.786944,
            )
            scenario = lin
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10

            assert np.all(scenario["S"] >= 0)
            assert np.all(scenario["I"] >= 0)

            cbr = np.random.uniform(5, 35, len(scenario))
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)

            R0 = 1.2
            infectious_duration_mean = 7.0
            beta = R0 / infectious_duration_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrate_map)

                infdist = dists.normal(loc=infectious_duration_mean, scale=2)
                pyramid = AliasedDistribution(np.full(89, 1_000))
                survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

                s = SIS.Susceptible(model)
                i = SIS.Infectious(model, infdist)
                tx = SIS.Transmission(model, infdist)
                births = BirthsByCBR(model, birthrate_map, pyramid)
                mortality = MortalityByEstimator(model, survival)
                model.components = [s, i, tx, births, mortality]

                model.validating = VALIDATING

            model.run(f"SIS Linear ({model.people.count:,}/{model.nodes.count:,})")

        # --- Validation checks ---
        assert np.all(model.nodes.S >= 0)
        assert np.all(model.nodes.I >= 0)
        assert np.allclose(
            (model.nodes.S + model.nodes.I)[0, :],
            scenario["population"],
            rtol=0.05,
        ), "Population mismatch at start"
        assert model.nodes.I[-1].sum() >= 0, "Infected counts must remain non-negative"

        if VERBOSE:
            print(model.people.describe("People"))
            print(model.nodes.describe("Nodes"))

        if PLOTTING:
            ibm = np.random.choice(len(base_maps))
            model.basemap_provider = base_maps[ibm]
            print(f"Using basemap: {model.basemap_provider.name}")
            model.plot()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--plot", action="store_true", help="Enable plotting")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("-m", type=int, default=5, help="Number of grid rows (M)")
    parser.add_argument("-n", type=int, default=5, help="Number of grid columns (N)")
    parser.add_argument("-p", type=int, default=10, help="Number of linear nodes (N)")
    parser.add_argument("--validating", action="store_true", help="Enable validating mode")
    parser.add_argument("-t", "--ticks", type=int, default=365, help="Number of days to simulate (nticks)")
    parser.add_argument("-g", "--grid", action="store_true", help="Run grid test")
    parser.add_argument("-l", "--linear", action="store_true", help="Run linear test")
    parser.add_argument("-c", "--constant", action="store_true", help="Run constant population test")
    parser.add_argument("unittest", nargs="*")
    args = parser.parse_args()

    PLOTTING = args.plot
    VERBOSE = args.verbose
    VALIDATING = args.validating
    NTICKS = args.ticks
    EM, EN, PEE = args.m, args.n, args.p

    print(f"Using arguments {args=}")

    if not (args.grid or args.linear or args.constant):
        sys.argv[1:] = args.unittest
        unittest.main(exit=False)
    else:
        tc = Default()
        if args.grid:
            tc.test_grid()
        if args.linear:
            tc.test_linear()

    ts.freeze()
    print("\nTiming Summary:")
    print("-" * 30)
    print(ts.to_string(scale="ms"))
    with Path("timing_data.json").open("w") as f:
        json.dump(ts.to_dict(scale="ms"), f, indent=4)
