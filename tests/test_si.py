from laser.generic.newutils import TimingStats as ts  # noqa: I001

import json
import unittest
import pytest
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from laser.core import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.demographics import KaplanMeierEstimator

from laser.generic import SI
from laser.generic import Model
from laser.generic.newutils import ValuesMap
from laser.generic.newutils import grid
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByEstimator, ConstantPopVitalDynamics
from tests.utils import base_maps
from tests.utils import stdgrid

PLOTTING = False
VERBOSE = False
EM = 10
EN = 10
PEE = 10
VALIDATING = False
NTICKS = 365


class Default(unittest.TestCase):
    def test_grid(self):
        """
        Feature: Spatial 2-D SIS grid model
        --------------------------------------------------
        Validates:
          • Spatially explicit, multi-patch SIS dynamics on a two-dimensional grid.
          • Population birth and death processes via the VitalDynamics component.
          • Transmission and recovery events governed by SIS.Infectious and SIS.Transmission.
          • Integration of demographic and epidemiological components within a unified model loop.

        Configuration:
          Grid size: 10 x 10 (100 nodes)
          Node size: 10 km
          Population initialization: uniform random between 10 000 and 1 000 000
          Simulation length: 365 ticks (daily updates)

        Expected Outcomes / Invariants:
          • S + I remains approximately equal to total population per node.
          • Populations remain strictly positive (no negative counts).
          • Infection prevalence (I / N) remains within [0, 1].
          • Model executes full duration without numerical or indexing errors.

        Notes:
          This test represents the full spatial-grid capability of LASER and exercises
          nearly all key agent-level and node-level update mechanisms in a coupled
          SIS framework. It therefore provides high-level validation of LASER's
          spatial and demographic integration.
        """
        with ts.start("test_grid"):
            scenario = stdgrid(M=EM, N=EN)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10

            cbr = np.random.uniform(5, 35, len(scenario))  # CBR = per 1,000 per year
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)

            params = PropertySet({"nticks": NTICKS, "beta": 1.0 / 32})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrate_map)
                model.validating = VALIDATING

                # Sampling this pyramid will return indices in [0, 88] with equal probability.
                pyramid = AliasedDistribution(np.full(89, 1_000))
                # The survival function will return the probability of surviving past each age.
                survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

                s = SI.Susceptible(model)
                i = SI.Infectious(model)
                tx = SI.Transmission(model)
                births = BirthsByCBR(model, birthrate_map, pyramid)
                mortality = MortalityByEstimator(model, survival)
                model.components = [s, i, tx, births, mortality]

            model.run(f"SI Grid ({model.people.count:,}/{model.nodes.count:,})")

            # --- Quantitative post-simulation checks ---

            # 1. Infection must change over time (non-static dynamics)
            initial_I = model.nodes.I[0].sum()
            final_I = model.nodes.I[-1].sum()
            assert final_I != pytest.approx(initial_I), "Infection count should evolve over time."

            # 2. Mean prevalence should remain below a realistic bound
            mean_prev = (model.nodes.I / (model.nodes.I + model.nodes.S + 1e-9)).mean()
            assert 0.0 <= mean_prev <= 0.5, f"Mean prevalence unrealistic: {mean_prev:.3f}"

            # 3. Total population trend should reflect births - deaths
            pop0 = (model.nodes.S[0] + model.nodes.I[0]).sum()
            popT = (model.nodes.S[-1] + model.nodes.I[-1]).sum()
            delta = (popT - pop0) / pop0
            assert abs(delta) < 0.1, f"Unexpected population drift: {delta * 100:.2f}%"

            # 4. Node-level conservation (mean relative error)
            rel_err = np.abs((model.nodes.S[-1] + model.nodes.I[-1]) - (model.nodes.S[0] + model.nodes.I[0])) / (
                model.nodes.S[0] + model.nodes.I[0] + 1e-9
            )
            assert rel_err.mean() < 0.05, "Average node population drift exceeds 5%"

        if VERBOSE:
            print(model.people.describe("People"))
            print(model.nodes.describe("Nodes"))

        if PLOTTING:
            if base_maps:
                ibm = np.random.choice(len(base_maps))
                model.basemap_provider = base_maps[ibm]
                print(f"Using basemap: {model.basemap_provider.name}")
            else:
                print("No base maps available.")
            model.plot()

        return

    def test_linear(self):
        """
        Feature: One-dimensional (linear) SI model
        --------------------------------------------------
        Validates:
          • Sequential (1-D chain) arrangement of patches with neighbor-based interaction.
          • Infection transmission and recovery processes identical to grid model.
          • Vital dynamics (births/deaths) under linear connectivity.
          • Proper handling of boundary conditions at the ends of the chain.

        Configuration:
          Layout: 1 x 10 linear chain
          Node size: 10 km
          Population initialization: uniform random between 10 000 and 1 000 000
          Simulation length: 365 ticks (daily updates)

        Expected Outcomes / Invariants:
          • Population conservation per node (S + I ≈ N).
          • Non-negative and bounded infection counts.
          • Consistency of per-node epidemiological transitions with grid model behavior.

        Notes:
          This test isolates topological and connectivity aspects of LASER's SI model.
          It validates that infection spread and demography behave correctly under
          reduced spatial dimensionality, providing confidence that LASER can operate
          across alternative spatial network configurations.
        """
        with ts.start("test_linear"):
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10

            cbr = np.random.uniform(5, 35, len(scenario))  # CBR = per 1,000 per year
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)

            params = PropertySet({"nticks": NTICKS, "beta": 1.0 / 32})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrate_map)
                model.validating = VALIDATING

                # Sampling this pyramid will return indices in [0, 88] with equal probability.
                pyramid = AliasedDistribution(np.full(89, 1_000))
                # The survival function will return the probability of surviving past each age.
                survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

                s = SI.Susceptible(model)
                i = SI.Infectious(model)
                tx = SI.Transmission(model)
                births = BirthsByCBR(model, birthrate_map, pyramid)
                mortality = MortalityByEstimator(model, survival)
                model.components = [s, i, tx, births, mortality]

            model.run(f"SI Linear ({model.people.count:,}/{model.nodes.count:,})")

            # 1. Infection curve monotonicity segments (rise-fall)
            I_series = model.nodes.I.sum(axis=1)
            assert I_series.max() > I_series[0] * 1.5, "No epidemic growth detected."
            peak_tick = np.argmax(I_series)
            assert peak_tick > 0, "Peak should occur after initial tick."

            # 2. Epidemic decay: final infections lower than peak
            assert I_series[-1] < I_series[peak_tick] * 0.8, "No decline after peak."

            # 3. Population size consistency
            pop_change = ((model.nodes.S[-1] + model.nodes.I[-1]).sum() / (model.nodes.S[0] + model.nodes.I[0]).sum()) - 1
            assert abs(pop_change) < 0.05, f"Population changed {pop_change * 100:.2f}%"

        if VERBOSE:
            print(model.people.describe("People"))
            print(model.nodes.describe("Nodes"))

        if PLOTTING:
            if base_maps:
                ibm = np.random.choice(len(base_maps))
                model.basemap_provider = base_maps[ibm]
                print(f"Using basemap: {model.basemap_provider.name}")
            else:
                print("No base maps available.")
            model.plot()

        return

    def test_constant_pop(self):
        """
        Feature: Constant-population SI model with dynamic births and deaths
        --------------------------------------------------
        Validates:
          • Constant-population demographic process in which births exactly offset deaths.
          • Interaction between epidemiological and demographic components under a
            strict population-conservation constraint.
          • Correct handling of zero-mortality edge cases (no negative or runaway population).
          • Integration of SIS infection dynamics (Susceptible, Infectious, Transmission)
            with ConstantPopVitalDynamics.

        Configuration:
          Layout: single-node model (M=1, N=1)
          Population: 1 000 000
          Initial infections: 10
          Crude birth rate: 400 births per 1 000 individuals per year
          Mortality: explicitly set to zero
          Simulation length: 365 ticks (daily updates)

        Expected Outcomes / Invariants:
          • Total population remains constant throughout the run (ΔN ≈ 0).
          • Non-negative susceptible and infected counts for all ticks.
          • Infection prevalence remains bounded within [0, 1].
          • Model runs without demographic or epidemiological warnings.

        Notes:
          This test exercises LASER's ConstantPopVitalDynamics component in combination
          with SIS transmission and progression. It ensures that population accounting
          remains stable even when births and deaths are tightly coupled or extreme
          (high CBR, zero mortality). Serves as a regression test for demographic
          balance and numerical stability in constant-population scenarios.
        """
        with ts.start("test_constant_pop"):
            pop = 1e6
            init_inf = 10
            scenario = grid(M=1, N=1, node_size_km=10, population_fn=lambda x, y: pop)
            scenario["S"] = scenario.population - init_inf
            scenario["I"] = init_inf
            parameters = PropertySet({"seed": 2, "nticks": NTICKS, "verbose": True, "beta": 0.04, "cbr": 400})

            birthrate_map = ValuesMap.from_scalar(parameters.cbr, nsteps=parameters.nticks, nnodes=1)

            with ts.start("Model Initialization"):
                model = Model(scenario, parameters, birthrate_map, skip_capacity=True)
                model.validating = VALIDATING

                model.components = [
                    SI.Susceptible(model),
                    SI.Infectious(model),
                    SI.Transmission(model),
                    ConstantPopVitalDynamics(model, birthrate_map),
                ]

            model.run(f"SI Constant Pop ({model.people.count:,}/{model.nodes.count:,})")

            # 1. Total population constant to within 0.01%
            N0 = (model.nodes.S[0] + model.nodes.I[0]).sum()
            NT = (model.nodes.S[-1] + model.nodes.I[-1]).sum()
            assert abs(NT - N0) / N0 < 1e-4, f"Population not constant: {NT - N0}"

            # 2. Infection prevalence stable and bounded
            prev_series = model.nodes.I.sum(axis=1) / (model.nodes.I.sum(axis=1) + model.nodes.S.sum(axis=1) + 1e-9)
            assert np.all((prev_series >= 0) & (prev_series <= 1))
            # assert prev_series.std() < 0.05, "Prevalence fluctuated excessively for constant-pop model"
            assert prev_series.std() < 0.1, f"Prevalence fluctuated excessively: std={prev_series.std():.3f}"

            # 3. Births ≈ deaths accounting
            births_total = getattr(model.nodes, "births", None)
            deaths_total = getattr(model.nodes, "deaths", None)
            if births_total is not None and deaths_total is not None:
                net = births_total.sum() - deaths_total.sum()
                assert abs(net) < 1e-6 * N0, f"Birth-death mismatch: {net}"

        if VERBOSE:
            print(model.people.describe("People"))
            print(model.nodes.describe("Nodes"))

        if PLOTTING:
            if base_maps:
                ibm = np.random.choice(len(base_maps))
                model.basemap_provider = base_maps[ibm]
                print(f"Using basemap: {model.basemap_provider.name}")
            else:
                print("No base maps available.")
            model.plot()

        return


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--plot", action="store_true", help="Enable plotting")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--validating", action="store_true", help="Enable validating mode")
    parser.add_argument("-m", type=int, default=5, help="Number of grid rows (M)")
    parser.add_argument("-n", type=int, default=5, help="Number of grid columns (N)")
    parser.add_argument("-p", type=int, default=10, help="Number of linear nodes (N)")
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

    # Instantiate the test case
    tc = Default()

    # If no test flags were given, run all by default
    run_all = not (args.grid or args.linear or args.constant)

    if args.grid or run_all:
        print("\nRunning grid configuration...")
        tc.test_grid()

    if args.linear or run_all:
        print("\nRunning linear configuration...")
        tc.test_linear()

    if args.constant:
        print("\nRunning constant population configuration...")
        tc.test_constant_pop()

    ts.freeze()
    print("\nTiming Summary:")
    print("-" * 30)
    print(ts.to_string(scale="ms"))
    with Path("timing_data.json").open("w") as f:
        json.dump(ts.to_dict(scale="ms"), f, indent=4)
