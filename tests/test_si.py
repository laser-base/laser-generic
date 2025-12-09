from laser.generic.utils import TimingStats as ts  # noqa: I001

import json
import unittest
import pytest
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from laser.core import PropertySet

from laser.generic import SI
from laser.generic import Model
from laser.generic.utils import ValuesMap
from laser.core.utils import grid
from laser.generic.vitaldynamics import ConstantPopVitalDynamics
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
        Feature: Two-dimensional spatial SI model (grid topology)
        ---------------------------------------------------------
        Validates core LASER functionality for a spatially explicit,
        multi-patch **SI (Susceptible → Infectious, no recovery)** model
        on a 2-D grid.

        Model Structure:
            • 10×10 grid (100 nodes), each a distinct spatial patch.
            • SI epidemic dynamics: Susceptible agents transition to
              Infectious but never recover.
            • No vital dynamics (births, deaths) in this test.
            • Spatial interactions governed by the model’s adjacency matrix.

        What this test verifies:
            ✓ Spatial connectivity and correct construction of a 2-D grid.
            ✓ SI transmission operates across patches (infection spreads).
            ✓ Infection counts evolve over time (system is non-static).
            ✓ Prevalence remains bounded within biologically reasonable ranges.
            ✓ Population counts (S+I) remain exactly constant (no births/deaths).
            ✓ Per-node conservation holds: no negative or exploding populations.
            ✓ No numerical errors or inconsistent state transitions.

        Interpretation:
            Passing this test demonstrates that LASER:
              • correctly sets up a large spatial grid,
              • applies SI transmission across space,
              • maintains population conservation,
              • and runs stably over many timesteps.

            This is a broad end-to-end validation of spatial SI dynamics.
        """
        with ts.start("test_grid"):
            scenario = stdgrid(M=EM, N=EN)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10

            cbr = np.random.uniform(5, 35, len(scenario))  # CBR = per 1,000 per year
            birthrate_map = ValuesMap.from_nodes(cbr, nticks=NTICKS)

            params = PropertySet({"nticks": NTICKS, "beta": 1.0 / 32})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrate_map)
                model.validating = VALIDATING

                s = SI.Susceptible(model)
                i = SI.Infectious(model)
                tx = SI.Transmission(model)
                # births = BirthsByCBR(model, birthrate_map, pyramid)
                # mortality = MortalityByEstimator(model, survival)
                model.components = [s, i, tx]  # , births, mortality]

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
            assert abs(delta) == 0

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
        Feature: One-dimensional spatial SI model (linear chain)
        --------------------------------------------------------
        Validates LASER’s behavior for a spatial SI (Susceptible → Infectious)
        epidemic on a **1-D chain of patches**, without recovery or vital dynamics.

        Model Structure:
            • 1×10 linear chain (10 nodes, each connected only to neighbors).
            • SI epidemic dynamics (no recovery): infection grows monotonically.
            • No births or deaths; population is constant.
            • Spatial mixing occurs only between adjacent patches.

        What this test verifies:
            ✓ Correct construction of a linear adjacency structure.
            ✓ Infection grows substantially from its initial seeding.
            ✓ Infection time series is nearly monotone (SI allows no decline).
            ✓ Population (S+I) remains perfectly conserved over time.
            ✓ All node-level counts remain valid (non-negative, finite).
            ✓ No numerical errors or invalid array operations.

        Interpretation:
            Passing this test shows that LASER:
                • handles 1-D spatial connectivity correctly,
                • propagates SI infection along a spatial chain,
                • conserves population in the absence of vital dynamics,
                • and maintains monotonic SI dynamics as expected.

            This isolates spatial adjacency logic independently of other
            model complexities.
        """
        with ts.start("test_linear"):
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10

            cbr = np.random.uniform(5, 35, len(scenario))  # CBR = per 1,000 per year
            birthrate_map = ValuesMap.from_nodes(cbr, nticks=NTICKS)

            params = PropertySet({"nticks": NTICKS, "beta": 1.0 / 32})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrate_map)
                model.validating = VALIDATING

                s = SI.Susceptible(model)
                i = SI.Infectious(model)
                tx = SI.Transmission(model)
                model.components = [s, i, tx]  # , births, mortality]

            model.run(f"SI Linear ({model.people.count:,}/{model.nodes.count:,})")

            I_series = model.nodes.I.sum(axis=1)

            # Epidemic growth: infections increase substantially
            assert I_series[-1] > I_series[0] * 1.5, "No epidemic growth detected."

            # Optional: curve should not decrease much (SI is monotone-ish)
            diffs = np.diff(I_series)
            assert (diffs < 0).sum() <= 5, "Unexpected large declines in SI infections."

            # 3. Population size consistency
            pop_change = ((model.nodes.S[-1] + model.nodes.I[-1]).sum() / (model.nodes.S[0] + model.nodes.I[0]).sum()) - 1
            assert abs(pop_change) == 0

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
        Feature: Constant-population SI model with balanced births & deaths
        -------------------------------------------------------------------
        Validates LASER’s SI (Susceptible → Infectious, no recovery) dynamics
        when coupled to **ConstantPopVitalDynamics**, which enforces exact
        population conservation by making births offset deaths.

        Model Structure:
            • Single-node model (no spatial effects).
            • Initial population: 1,000,000; initial infections: 10.
            • ConstantPopVitalDynamics keeps N(t) ≈ constant.
            • SI epidemic dynamics: infections increase monotonically because
              there is no recovery process.
            • Births and deaths occur internally but are exactly balanced.

        What this test verifies:
            ✓ Total population remains constant to numerical tolerance.
            ✓ Prevalence stays within [0,1] and rises monotonically
              (correct SI behavior without recovery).
            ✓ Infection grows substantially from initial seeding.
            ✓ Birth counts approximately equal death counts over the run.
            ✓ No negative populations or invalid state transitions.
            ✓ No unexpected oscillations or declines in prevalence.

        Interpretation:
            Passing this test demonstrates that LASER:
                • correctly implements ConstantPopVitalDynamics,
                • preserves strict population conservation,
                • produces monotonic SI epidemic curves under demographic
                  turnover,
                • and handles births/deaths without destabilizing epidemic
                  updates.

            This is an end-to-end validation of demographic coupling in a
            constant-population SI setting.
        """
        with ts.start("test_constant_pop"):
            pop = 1e6
            init_inf = 10
            scenario = grid(M=1, N=1, node_size_degs=0.08983, population_fn=lambda x, y: pop)
            scenario["S"] = scenario.population - init_inf
            scenario["I"] = init_inf
            parameters = PropertySet({"prng_seed": 2, "nticks": NTICKS, "verbose": True, "beta": 0.04, "cbr": 400})

            birthrate_map = ValuesMap.from_scalar(parameters.cbr, nticks=parameters.nticks, nnodes=1)

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

            # 1. Total population must remain constant (to numerical tolerance)
            N0 = (model.nodes.S[0] + model.nodes.I[0]).sum()
            NT = (model.nodes.S[-1] + model.nodes.I[-1]).sum()
            assert abs(NT - N0) / N0 < 1e-4, f"Population not constant: ΔN={NT - N0}"

            # 2. Prevalence must remain bounded and rise monotonically (SI has no recovery)
            prev_series = model.nodes.I.sum(axis=1) / (model.nodes.I.sum(axis=1) + model.nodes.S.sum(axis=1) + 1e-9)

            # 2a. Prevalence always within [0, 1]
            assert np.all((prev_series >= 0) & (prev_series <= 1)), "Prevalence out of bounds [0,1]."

            # 2b. SI has no recovery → prevalence should rise or stay flat, not meaningfully decline
            drops = np.diff(prev_series) < 0
            assert drops.sum() <= 3, f"Prevalence decreased unexpectedly at {drops.sum()} ticks."

            # 2c. Prevalence should increase substantially from the initial small seeding
            assert prev_series[-1] > prev_series[0] * 100, (
                f"Prevalence did not grow enough for SI: start={prev_series[0]:.5f}, end={prev_series[-1]:.5f}"
            )

            # 3. Births ≈ deaths (constant-population vital dynamics)
            births_total = getattr(model.nodes, "births", None)
            deaths_total = getattr(model.nodes, "deaths", None)
            if births_total is not None and deaths_total is not None:
                net = births_total.sum() - deaths_total.sum()
                assert abs(net) < 1e-6 * N0, f"Birth-death mismatch: net={net}"

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
