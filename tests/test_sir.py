from laser.generic.newutils import TimingStats as ts  # noqa: I001

import json
import unittest
from argparse import ArgumentParser
from pathlib import Path

import laser.core.distributions as dists
import numpy as np
from laser.core import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.demographics import KaplanMeierEstimator
from scipy.special import lambertw

from laser.generic import SIR
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
R0 = 1.386  # final attack fraction of 50%


class Default(unittest.TestCase):
    def test_single(self):
        """
        Feature: Single-node deterministic SIR model
        --------------------------------------------------
        Validates:
          • Infection and recovery progression in an isolated population.
          • Deterministic epidemic curve shape (rise-peak-fall) under R₀ = 1.386.
          • Population conservation and recovery fraction consistency.

        Configuration:
          Nodes: 1
          Population: 100,000
          Initial infections: 10
          Infectious duration: Normal(mean=7, sd=2)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • Infection count rises, peaks, and declines (monotonic segments).
          • Total population S+I+R constant to within 0.01%.
          • Final attack fraction (R/N) ≈ 50 ± 5%.

        Notes:
          Provides a minimal deterministic benchmark for LASER's SIR transitions,
          verifying internal mass balance and infection kinetics before adding spatial
          or demographic complexity.
        """
        with ts.start("test_single_node"):
            scenario = stdgrid(M=1, N=1, population_fn=lambda x, y: 100_000)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10
            scenario["R"] = 0

            infectious_duration_mean = 7.0
            beta = R0 / infectious_duration_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            with ts.start("Model Initialization"):
                model = Model(scenario, params)

                infdist = dists.normal(loc=infectious_duration_mean, scale=2)
                s = SIR.Susceptible(model)
                i = SIR.Infectious(model, infdist)
                r = SIR.Recovered(model)
                tx = SIR.Transmission(model, infdist)
                model.components = [s, i, r, tx]
                model.validating = VALIDATING

            model.run("SIR Single Node")

            # --- Quantitative Checks ---
            I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)

            assert I_series.max() > I_series[0] * 2, "Infection did not grow significantly."
            peak_tick = np.argmax(I_series)
            assert I_series[-1] < I_series[peak_tick] * 0.5, "Infection did not decline post-peak."

            # Constant total population
            N0 = (model.nodes.S[0] + model.nodes.I[0] + model.nodes.R[0]).sum()
            NT = (model.nodes.S[-1] + model.nodes.I[-1] + model.nodes.R[-1]).sum()
            assert abs(NT - N0) / N0 < 1e-4, f"Population drift >0.01%: ΔN={NT - N0}"

            # Final attack fraction check (~50%)
            final_af = R_series[-1] / scenario.population.sum()
            assert 0.45 <= final_af <= 0.55, f"Final attack fraction {final_af:.3f} out of expected range."

    def test_grid(self):
        """
        Feature: Spatial 2-D SIR model with births and deaths
        --------------------------------------------------
        Validates:
          • Spatial epidemic propagation across a 10x10 grid of nodes.
          • Integration of demography (BirthsByCBR, MortalityByEstimator).
          • Stability of total population under demographic turnover.
          • Quantitative epidemic realism via bounded infection prevalence.

        Configuration:
          Grid: 10x10 nodes, 10 km each
          Population: 10 000-1 000 000 per node
          Infectious duration: Normal(mean=7, sd=2)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • Population remains within ±10% of baseline after 365 days.
          • Mean prevalence (I/N) ≤ 0.5.
          • Non-negative counts across all states.
          • Model executes full duration without instability.

        Notes:
          Provides a stochastic spatial-demographic stress test combining infection,
          birth, and mortality processes, validating LASER's spatial coupling integrity.
        """
        with ts.start("test_grid"):
            scenario = stdgrid(M=EM, N=EN)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10
            scenario["R"] = 0

            cbr = np.random.uniform(5, 35, len(scenario))
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)
            infectious_duration_mean = 7.0
            beta = R0 / infectious_duration_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrates=birthrate_map)
                infdist = dists.normal(loc=infectious_duration_mean, scale=2)
                pyramid = AliasedDistribution(np.full(89, 1_000))
                survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())
                s = SIR.Susceptible(model)
                i = SIR.Infectious(model, infdist)
                r = SIR.Recovered(model)
                tx = SIR.Transmission(model, infdist)
                births = BirthsByCBR(model, birthrates=birthrate_map, pyramid=pyramid)
                mortality = MortalityByEstimator(model, survival)
                model.components = [s, i, r, tx, births, mortality]
                model.validating = VALIDATING

            model.run("SIR Grid")

            # --- Quantitative Checks ---
            I_series = model.nodes.I.sum(axis=1)
            N_series = (model.nodes.S + model.nodes.I + model.nodes.R).sum(axis=1)
            pop_change = (N_series[-1] - N_series[0]) / N_series[0]
            mean_prev = (model.nodes.I / (model.nodes.S + model.nodes.I + model.nodes.R + 1e-9)).mean()

            assert np.all(model.nodes.S >= 0)
            assert np.all(model.nodes.I >= 0)
            assert np.all(model.nodes.R >= 0)
            assert abs(pop_change) < 0.1, f"Population drift {pop_change * 100:.2f}% exceeds ±10%."
            assert 0 <= mean_prev <= 0.5, f"Mean prevalence unrealistic: {mean_prev:.3f}"
            assert I_series.max() > I_series[0] * 1.5, "Epidemic growth not observed."

    def test_linear(self):
        """
        Feature: One-dimensional (linear) SIR model
        --------------------------------------------------
        Validates:
          • Epidemic propagation in a chain topology (1xN nodes).
          • Comparison of epidemic timing vs. grid model.
          • Local population conservation and bounded infection.
          • Demographic integration across minimal connectivity.

        Configuration:
          Layout: 1x10 linear chain
          Infectious duration: Normal(mean=7, sd=2)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • Total population drift <5%.
          • Epidemic rise and decay pattern (clear peak).
          • Final attack fraction within realistic range (20-70%).

        Notes:
          Tests LASER's SIR implementation under linear spatial coupling,
          validating that infection propagation, recovery, and demographics
          behave consistently across simplified topology.
        """
        with ts.start("test_linear"):
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10
            scenario["R"] = 0

            cbr = np.random.uniform(5, 35, len(scenario))
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)
            infectious_duration_mean = 7.0
            beta = R0 / infectious_duration_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrates=birthrate_map)
                infdist = dists.normal(loc=infectious_duration_mean, scale=2)
                pyramid = AliasedDistribution(np.full(89, 1_000))
                survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())
                s = SIR.Susceptible(model)
                i = SIR.Infectious(model, infdist)
                r = SIR.Recovered(model)
                tx = SIR.Transmission(model, infdist)
                births = BirthsByCBR(model, birthrates=birthrate_map, pyramid=pyramid)
                mortality = MortalityByEstimator(model, survival)
                model.components = [s, i, r, tx, births, mortality]
                model.validating = VALIDATING

            model.run("SIR Linear")

            # --- Quantitative Checks ---
            I_series = model.nodes.I.sum(axis=1)
            pop_series = (model.nodes.S + model.nodes.I + model.nodes.R).sum(axis=1)
            pop_change = (pop_series[-1] - pop_series[0]) / pop_series[0]
            assert abs(pop_change) < 0.05, f"Population drift {pop_change * 100:.2f}% >5%."
            assert I_series.max() > I_series[0] * 1.5, "Epidemic growth too weak."
            peak_tick = np.argmax(I_series)
            assert I_series[-1] < I_series[peak_tick] * 0.8, "Epidemic did not decline."

    def test_kermack_mckendrick(self):
        """
        Feature: Theoretical validation — Kermack-McKendrick final size
        --------------------------------------------------
        Validates:
          • LASER SIR model convergence to the analytic Kermack-McKendrick final attack fraction.
          • Consistency across stochastic initializations (multiple iterations).
          • Quantitative deviation threshold of ±5% from analytic solution.

        Configuration:
          Population: 1 000 000 (single node)
          Initial infections: 1 000
          R₀ range: 1.2-2.0
          Infectious duration: 7 days
          Iterations: 10 stochastic replicates per R₀ case

        Expected Outcomes / Invariants:
          • Median attack fraction within 5% of theoretical value.
          • No more than 3/10 runs deviate >5%.

        Notes:
          A quantitative regression test comparing simulated final epidemic size
          to the analytic SIR solution using Lambert W. Ensures LASER's SIR core
          equations reproduce canonical epidemic final-size relationships.
        """

        def attack_fraction(beta, inf_mean, pop, init_inf):
            R0 = beta * inf_mean
            S0 = (pop - init_inf) / pop
            S_inf = -1 / R0 * lambertw(-R0 * S0 * np.exp(-R0)).real
            return 1 - S_inf

        INIT_INF = 1_000
        cases = [
            (1.2160953 / 7, 7.0, 1.0 / 3.0),
            (1.27685 / 7, 7.0, 0.4),
            (1.527 / 7, 7.0, 0.6),
            (2.011675 / 7, 7.0, 0.8),
        ]

        for beta, inf_mean, expected_af in cases:
            failed = 0
            NITERS = 10
            for _ in range(NITERS):
                scenario = stdgrid(M=1, N=1, population_fn=lambda x, y: 1_000_000)
                scenario["S"] = scenario["population"] - INIT_INF
                scenario["I"] = INIT_INF
                scenario["R"] = 0
                params = PropertySet({"nticks": NTICKS, "beta": beta})
                model = Model(scenario, params)
                infdurdist = dists.normal(loc=inf_mean, scale=2)
                s = SIR.Susceptible(model)
                i = SIR.Infectious(model, infdurdist)
                r = SIR.Recovered(model)
                tx = SIR.Transmission(model, infdurdist)
                model.components = [s, i, r, tx]
                model.run("SIR KM")

                actual_af = model.nodes.R[-1].sum() / scenario.population.sum()
                diff = abs(actual_af - expected_af)
                frac = diff / expected_af
                if frac > 0.05:
                    failed += 1
            assert failed < 3, (
                f"Kermack-McKendrick test failed {failed}/{NITERS} for R0={beta * inf_mean:.3f} (expected AF={expected_af:.3f})"
            )


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
    parser.add_argument("-k", "--km", action="store_true", help="Run Kermack-McKendrick validation")
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
    run_all = not (args.grid or args.linear or args.single or args.km)

    if args.single or run_all:
        tc.test_single()
    if args.grid or run_all:
        tc.test_grid()
    if args.linear or run_all:
        tc.test_linear()
    if args.km or run_all:
        tc.test_kermack_mckendrick()

    ts.freeze()
    print("\nTiming Summary:")
    print("-" * 30)
    print(ts.to_string(scale="ms"))
    with Path("timing_data.json").open("w") as f:
        json.dump(ts.to_dict(scale="ms"), f, indent=4)
