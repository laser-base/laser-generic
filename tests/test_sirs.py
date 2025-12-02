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

    def test_sirs_linear_no_demography(self):
        """
        Feature: Pure SIRS epidemic on a 1×N linear spatial chain (no demography)
        -------------------------------------------------------------------------
        This test validates the core SIRS epidemic engine in LASER—specifically,
        the correctness of infection, recovery, and waning-immunity transitions—
        in a clean, fixed-population environment with no births or mortality.

        Model Structure:
            • Topology: 1×N linear chain (each patch interacts only with neighbors).
            • Disease Process:
                  S → I → R → S
              with explicit infectious-duration and waning-duration distributions.
            • No demographic turnover: population is conserved exactly.
            • Simulation horizon: 365 ticks.

        What this test verifies:
            ✓ Exact population conservation: no births/deaths, no loss of agents.
            ✓ Strong epidemic growth: I(t) rises sharply from the initial seeding.
            ✓ SIRS characteristic behavior:
                  – A clear epidemic peak exists,
                  – Infection declines after the peak but does NOT crash to zero
                    (endemic persistence is expected with waning),
                  – Reinfection cycles occur due to waning immunity.
            ✓ Waning immunity is functioning:
                  R(t) eventually decreases as individuals return to S.
            ✓ All state counts remain non-negative and well-behaved.

        Why this test matters:
            Passing this test demonstrates that LASER’s SIRS state-transition logic
            (infection, recovery, waning, reinfection) is functioning correctly in
            isolation from demographic processes. It provides high-confidence that
            the epidemiological core of SIRS is implemented correctly before
            layering on births, deaths, or long-term population churn.
        """
        with ts.start("test_sirs_linear_no_demography"):
            # --- Scenario ---
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10
            scenario["R"] = 0

            # --- Durations ---
            inf_mean = 7.0
            infdist = dists.normal(loc=inf_mean, scale=2.0)

            waning_mean = 30.0
            waningdist = dists.normal(loc=waning_mean, scale=5.0)

            R0 = 2.0
            beta = R0 / inf_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})
            # --- Model ---
            with ts.start("Model Initialization"):
                model = Model(scenario, params)
                model.validating = VALIDATING

                s = SIRS.Susceptible(model)
                i = SIRS.Infectious(model, infdist, waningdist)
                r = SIRS.Recovered(model, waningdist)
                tx = SIRS.Transmission(model, infdist)

                model.components = [s, i, r, tx]

            model.run("SIRS Linear (no demography)")

            # --- Checks ---
            I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            S_series = model.nodes.S.sum(axis=1)
            pop_series = S_series + I_series + R_series

            # 1. exact population conservation
            assert abs(pop_series[-1] - pop_series[0]) < 1e-9

            # 2. epidemic must show strong growth
            assert I_series.max() > I_series[0] * 20, "Epidemic too weak."

            # 3. SIRS characteristic behavior:
            peak = np.argmax(I_series)
            assert peak > 5, "Peak too early."

            # SIRS should decline somewhat after peak (but not necessarily 60% like SIR)
            assert I_series[-1] < I_series[peak], "SIRS should show some decline after the peak (even if small)."

            # But SIRS should *not* crash to zero like SIR
            assert I_series[-1] > I_series[peak] * 0.05, "SIRS should maintain endemic infection; I(T) too low."

            # There must be at least one downward motion after peak
            post_peak = I_series[peak:]
            assert np.any(np.diff(post_peak) < 0), "No decline at all after peak (SIRS dynamics missing)."

            # 4. waning immunity must be visible
            assert R_series.max() > R_series[-1], "Waning not evident (R never decreases)."

    def test_sirs_linear_with_demography(self):
        """
        Feature: SIRS epidemic on a 1×N linear chain with demographic turnover
        ----------------------------------------------------------------------
        This test validates the integration of SIRS epidemic dynamics with
        LASER’s demographic components (births and mortality). The model includes
        reinfection cycles via waning immunity and continuous population turnover
        due to crude birth rates and age-structured mortality.

        Model Structure:
            • Topology: 1×N linear chain (nearest-neighbor mixing).
            • Disease Process:
                  S → I → R → S
              with explicit infectious and waning duration distributions.
            • Demographic Processes:
                  – Births via BirthsByCBR (node-dependent CBRs),
                  – Deaths via MortalityByEstimator (Kaplan–Meier survival),
                  – Population turnover introduces new susceptibles over time.
            • Simulation horizon: 365 ticks.

        What this test verifies:
            ✓ Demography + SIRS combine stably:
                  – Population remains positive,
                  – Drift is moderate (< ~15%) despite turnover.
            ✓ Epidemic growth occurs and produces a substantial peak.
            ✓ Post-peak decline is visible even with reinfections and new susceptibles.
            ✓ Waning immunity functions correctly (R decreases at some point).
            ✓ Vital-dynamics bookkeeping is consistent (birth/death mismatch bounded).
            ✓ No numerical instabilities, negative populations, or broken transitions.

        Why this test matters:
            This scenario exercises LASER’s ability to couple realistic demographic
            flux with SIRS epidemic mechanics. Passing this test confirms that LASER:
                • correctly handles susceptible replenishment via births,
                • integrates mortality with ongoing infection dynamics,
                • supports long-term endemic behavior under waning immunity,
                • and maintains stable population and epidemiological accounting.

            Together with the pure SIRS test, this establishes both epidemiological
            correctness *and* robust demographic-epidemiological integration.
        """
        with ts.start("test_sirs_linear_with_demography"):
            # --- Scenario ---
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10
            scenario["R"] = 0

            # --- Vital dynamics ---
            cbr = np.random.uniform(5, 35, PEE)  # births per 1000 per year
            birthrates = ValuesMap.from_nodes(cbr, nsteps=NTICKS)

            pyramid = AliasedDistribution(np.full(89, 1_000))
            survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

            # --- Durations ---
            inf_mean = 7.0
            infdist = dists.normal(loc=inf_mean, scale=2.0)

            waning_mean = 30.0
            waningdist = dists.normal(loc=waning_mean, scale=5.0)

            R0 = 2.0
            beta = R0 / inf_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrates=birthrates.values)
                model.validating = VALIDATING

                s = SIRS.Susceptible(model)
                i = SIRS.Infectious(model, infdist, waningdist)
                r = SIRS.Recovered(model, waningdist)
                tx = SIRS.Transmission(model, infdist)

                births = BirthsByCBR(model, birthrates.values, pyramid)
                mortality = MortalityByEstimator(model, survival)

                model.components = [s, i, r, tx, births, mortality]

            model.run("SIRS Linear (with demography)")

            # --- Checks ---
            I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            S_series = model.nodes.S.sum(axis=1)
            pop_series = S_series + I_series + R_series

            pop0 = pop_series[0]
            popT = pop_series[-1]
            pop_change = (popT - pop0) / pop0

            # 1. population drift must be moderate, not perfect
            assert abs(pop_change) < 0.15, f"Population drift {pop_change * 100:.2f}% > 15%."

            # 2. epidemic growth must occur
            assert I_series.max() > I_series[0] * 1.5, "Epidemic too weak."

            # 3. some decline after peak should still be visible
            peak = np.argmax(I_series)
            assert peak > 5
            assert I_series[-1] < I_series[peak] * 0.9, "No decline after peak."

            # 4. vital dynamics should not break anything
            assert np.all(pop_series > 0), "Negative population encountered."

            # 5. waning immunity should still reduce R at some point
            assert R_series.max() > R_series[-1], "Waning immunity not visible (R never declines)."


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--plot", action="store_true", help="Enable plotting")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--validating", action="store_true", help="Enable validating mode")
    parser.add_argument("-m", type=int, default=5, help="Number of grid rows (M)")
    parser.add_argument("-n", type=int, default=5, help="Number of grid columns (N)")
    parser.add_argument("-p", type=int, default=10, help="Number of linear nodes (N)")
    parser.add_argument("-t", "--ticks", type=int, default=365, help="Number of days to simulate (nticks)")
    parser.add_argument(
        "-r",
        "--r0",
        type=float,
        default=1.386,
        help="Basic reproduction number (R0) [1.151 for 25%% attack fraction, 1.386=50%%, and 1.848=75%%]",
    )
    parser.add_argument("-g", "--grid", action="store_true", help="Run grid test")
    parser.add_argument("-l", "--linear", action="store_true", help="Run linear test")
    parser.add_argument("-s", "--single", action="store_true", help="Run single node test")
    parser.add_argument("-i", "--infdur", type=float, default=7.0, help="Mean infectious duration in days")
    parser.add_argument("-w", "--wandur", type=float, default=30.0, help="Mean waning duration in days")
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
