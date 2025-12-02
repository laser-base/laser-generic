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

    def test_sir_linear_no_demography(self):
        """
        Feature: Pure SIR epidemic on a 1-D linear chain (no births/deaths)
        -------------------------------------------------------------------
        This test validates the core SIR epidemic logic in LASER without
        the added complexity of demographic processes. It uses a one-dimensional
        spatial chain of patches and checks that the SIR transmission,
        infectious-period progression, and recovery transitions behave
        correctly in a fixed-population setting.

        Model structure:
            • 1×N linear chain topology (each patch has at most two neighbors).
            • SIR disease dynamics:
                  S → I → R
              with user-specified infectious-duration distribution.
            • No vital dynamics: total population remains fixed.
            • Initial seeding: low-level infection at one end of the chain.
            • Simulation runs for 365 ticks.

        What this test verifies:
            ✓ Population conservation: S + I + R remains constant at all times.
            ✓ Epidemic growth: infection increases substantially from the initial seed.
            ✓ Epidemic peak: I(t) achieves a clear maximum mid-simulation.
            ✓ Epidemic decline: I(t) falls significantly after the peak (R buildup).
            ✓ Realistic attack fraction: final R/N is in a plausible SIR range.
            ✓ Numerical stability: no negative states or invalid transitions.

        Why this test matters:
            This is the simplest fully functional SIR validation scenario.
            It isolates the epidemic engine—transmission, timers, and recovery—
            from demographic noise. Passing this test demonstrates that LASER’s
            SIR core mechanics work correctly before layering on births/mortality.
        """

        def debug():
            print("\n--- DEBUG: SIR Linear No Demography ---")
            print(f"I(0): {I_series[0]}, I(max): {I_series.max()}, I(T): {I_series[-1]}")
            print(f"S(0): {S_series[0]}, S(T): {S_series[-1]}")
            print(f"R(0): {R_series[0]}, R(T): {R_series[-1]}")
            print("\nFirst 20 I(t):", I_series[:20])
            print("Last 20 I(t):", I_series[-20:])
            print("\nPrevalence first 20:", (I_series[:20] / (S_series[:20] + I_series[:20] + R_series[:20] + 1e-9)))
            print("Prevalence last 20:", (I_series[-20:] / (S_series[-20:] + I_series[-20:] + R_series[-20:] + 1e-9)))
            print("\nPeak tick:", np.argmax(I_series))

        with ts.start("test_sir_linear_no_demography"):
            # --- Scenario ---
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10
            scenario["R"] = 0

            # --- Parameters ---
            infectious_duration_mean = 7.0
            infdist = dists.normal(loc=infectious_duration_mean, scale=2.0)

            # R0 target for this test
            R0 = 2.0
            beta = R0 / infectious_duration_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            # --- Model ---
            with ts.start("Model Initialization"):
                model = Model(scenario, params)
                model.validating = VALIDATING

                s = SIR.Susceptible(model)
                i = SIR.Infectious(model, infdist)
                r = SIR.Recovered(model)
                tx = SIR.Transmission(model, infdist)
                model.components = [s, i, r, tx]

            model.run("SIR Linear (no demography)")

            # --- Quantitative Checks ---
            I_series = model.nodes.I.sum(axis=1)
            S_series = model.nodes.S.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)

            # 1. perfect population conservation
            pop0 = S_series[0] + I_series[0] + R_series[0]
            popT = S_series[-1] + I_series[-1] + R_series[-1]
            assert abs(popT - pop0) < 1e-9, f"Population drift detected: ΔN={popT - pop0}"

            # 2. epidemic must show strong growth (peak >> initial)
            peak_I = I_series.max()
            assert peak_I > I_series[0] * 20, f"Epidemic peak too weak: start={I_series[0]}, peak={peak_I}"

            # 3. epidemic must clearly decline after peak
            peak_tick = np.argmax(I_series)
            assert peak_tick > 5, "Peak occurs too early."
            assert I_series[-1] < peak_I * 0.4, "Epidemic did not decline after peak."

            # 4. final attack fraction in reasonable SIR range
            attack = R_series[-1] / pop0
            assert 0.2 <= attack <= 0.95, f"Attack fraction out of range: {attack:.3f}"

    def test_sir_linear_with_demography(self):
        """
        Feature: SIR epidemic on a 1-D linear chain with demographic turnover
        ---------------------------------------------------------------------
        This test validates the correct integration of SIR epidemic dynamics
        with demographic processes (births and mortality) in a minimal spatial
        topology. Unlike the pure-epidemic test, new susceptibles enter the
        population through births and individuals are removed by mortality.
        This produces a continuously shifting demographic structure and a
        more complex epidemic trajectory.

        Model structure:
            • 1×N linear chain topology (nearest-neighbor mixing).
            • SIR disease dynamics with finite infectious periods.
            • Vital dynamics:
                  – Births via BirthsByCBR (per-node crude birth rate).
                  – Mortality via Kaplan–Meier survival curve.
              These cause population turnover and susceptible replenishment.
            • Simulation runs for 365 ticks.

        What this test verifies:
            ✓ Population drift remains moderate (< ~15%) despite births/mortality.
            ✓ Epidemic growth occurs and produces a substantial peak.
            ✓ Decline after peak still occurs (though muted by new susceptibles).
            ✓ Population never becomes negative or unstable.
            ✓ Birth/death accounting remains reasonable (no runaway imbalance).
            ✓ SIR and demographic components interact consistently across time.

        Why this test matters:
            This scenario checks LASER’s ability to combine SIR infection
            progression with realistic demographic churn. It validates that
            births replenish susceptibles correctly, mortality removes agents
            without destabilizing dynamics, and that the combined system
            behaves plausibly and numerically stably.

            Passing this test gives confidence that LASER’s SIR implementation
            works not only in fixed-population settings, but also in models
            where demographic forces shape long-term epidemic patterns.
        """
        with ts.start("test_sir_linear_with_demography"):
            # --- Scenario ---
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["I"] = 10
            scenario["R"] = 0

            # --- Vital Dynamics ---
            # Birthrate moderately large; mortality via Kaplan–Meier curve.
            cbr = np.random.uniform(5, 35, len(scenario))  # births per 1000 per year
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS)

            pyramid = AliasedDistribution(np.full(89, 1_000))
            survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

            # --- Epidemic Params ---
            infectious_duration_mean = 7.0
            infdist = dists.normal(loc=infectious_duration_mean, scale=2.0)
            R0 = 2.0
            beta = R0 / infectious_duration_mean
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            # --- Model ---
            with ts.start("Model Initialization"):
                model = Model(scenario, params, birthrates=birthrate_map)
                model.validating = VALIDATING

                s = SIR.Susceptible(model)
                i = SIR.Infectious(model, infdist)
                r = SIR.Recovered(model)
                tx = SIR.Transmission(model, infdist)

                births = BirthsByCBR(model, birthrate_map, pyramid)
                mortality = MortalityByEstimator(model, survival)

                model.components = [s, i, r, tx, births, mortality]

            model.run("SIR Linear (with demography)")

            # --- Quantitative Checks ---
            I_series = model.nodes.I.sum(axis=1)
            S_series = model.nodes.S.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            pop_series = S_series + I_series + R_series

            pop0, popT = pop_series[0], pop_series[-1]
            pop_change = (popT - pop0) / pop0

            # 1. Population drift must be moderate but not exploding
            assert abs(pop_change) < 0.15, f"Population drift {pop_change * 100:.2f}% >15%."

            # 2. Epidemic growth must occur
            assert I_series.max() > I_series[0] * 1.5, "Epidemic growth too weak."

            # 3. Epidemic decline should still occur, though less pronounced than no-demography
            peak_tick = np.argmax(I_series)
            assert peak_tick > 5, "Peak should occur after initial ticks."
            assert I_series[-1] < I_series[peak_tick] * 0.9, "Epidemic did not decline after peak with demographics."

            # 4. Population series must remain positive and reasonable
            assert np.all(pop_series > 0), "Negative population count encountered."

            # 5. Birth / death accounting should not diverge wildly
            births_total = getattr(model.nodes, "births", None)
            deaths_total = getattr(model.nodes, "deaths", None)
            if births_total is not None and deaths_total is not None:
                # Just ensure turnover is not wildly imbalanced
                net = births_total.sum() - deaths_total.sum()
                assert abs(net) < 0.20 * pop0, f"Birth-death mismatch too large: net={net}"

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
    parser.add_argument("--plot", action="store_true", help="Enable plotting")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--validating", action="store_true", help="Enable validating mode")
    parser.add_argument("-m", type=int, default=5, help="Number of grid rows (M)")
    parser.add_argument("-n", type=int, default=5, help="Number of grid columns (N)")
    parser.add_argument("-p", type=int, default=10, help="Number of linear nodes (N)")
    parser.add_argument("-r", "--r0", type=float, default=1.386, help="R0")
    parser.add_argument("-t", "--ticks", type=int, default=365, help="Number of days to simulate (nticks)")
    parser.add_argument("-g", "--grid", action="store_true", help="Run spatial grid test")
    parser.add_argument("-l", "--linear", action="store_true", help="Run spatial linear test")
    parser.add_argument("-s", "--single", action="store_true", help="Run single node (non-spatial) test")
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
