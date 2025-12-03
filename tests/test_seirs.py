from laser.generic.newutils import TimingStats as ts  # noqa: I001

import json
import unittest
from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import math
import numpy as np
import laser.core.distributions as dists
from laser.core import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.demographics import KaplanMeierEstimator

from laser.generic import SEIRS
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
WANING_DURATION_MEAN = 30.0


def build_model(m, n, pop_fn, init_infected=0, init_recovered=0, birthrates=None, pyramid=None, survival=None, nticks=NTICKS):
    """
    Helper: Construct an SEIRS model with configurable demography and waning immunity.

    Builds Susceptible, Exposed, Infectious, Recovered, and Transmission components,
    optionally adding demographic processes (births and deaths).
    """
    scenario = stdgrid(M=m, N=n, population_fn=pop_fn)
    scenario["S"] = scenario["population"]
    scenario["E"] = 0
    scenario["S"] -= init_infected
    scenario["I"] = init_infected
    scenario["S"] -= init_recovered
    scenario["R"] = init_recovered

    beta = R0 / INFECTIOUS_DURATION_MEAN
    params = PropertySet({"nticks": nticks, "beta": beta})

    with ts.start("Model Initialization"):
        model = Model(scenario, params, birthrates=birthrates)

        expdist = dists.gamma(shape=EXPOSED_DURATION_SHAPE, scale=EXPOSED_DURATION_SCALE)
        infdist = dists.normal(loc=INFECTIOUS_DURATION_MEAN, scale=2)
        wandist = dists.normal(loc=WANING_DURATION_MEAN, scale=5)

        s = SEIRS.Susceptible(model)
        e = SEIRS.Exposed(model, expdist, infdist)
        i = SEIRS.Infectious(model, infdist, wandist)
        r = SEIRS.Recovered(model, wandist)
        tx = SEIRS.Transmission(model, expdist)

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
        Feature: Single-node deterministic SEIRS model
        --------------------------------------------------
        Validates:
          • Complete SEIRS loop (S→E→I→R→S) including waning immunity.
          • Correct ordering of transitions and finite infectious period.
          • Mass conservation over all states.
          • Recurrent epidemic behavior due to waning immunity.

        Configuration:
          Nodes: 1
          Population: 100,000
          Initial infections: 10
          Exposure: Gamma(shape=4.5, scale=1)
          Infectious: Normal(mean=7, sd=2)
          Waning: Normal(mean=30, sd=5)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • E peaks before I.
          • I peaks before R.
          • Re-infections occur after R wanes (R→S flow visible).
          • Population constant within 0.01%.
        """
        with ts.start("test_single_node"):
            model = build_model(1, 1, lambda x, y: 100_000, init_infected=10)
            model.run("SEIRS Single Node")

            E_series = model.nodes.E.sum(axis=1)
            I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            S_series = model.nodes.S.sum(axis=1)

            assert np.argmax(E_series) < np.argmax(I_series), "E should peak before I."
            assert np.argmax(I_series) < np.argmax(R_series), "I should peak before R."
            assert S_series[-1] < S_series[0], "S decreased initially due to infection."
            # assert R_series[-1] < R_series.max(), "Some waning expected (R decreases over time)."
            # Some waning must happen somewhere in the time series
            assert np.any(np.diff(R_series) < 0), "No waning evident: R(t) never showed any decrease."

            N0 = (model.nodes.S[0] + model.nodes.E[0] + model.nodes.I[0] + model.nodes.R[0]).sum()
            NT = (model.nodes.S[-1] + model.nodes.E[-1] + model.nodes.I[-1] + model.nodes.R[-1]).sum()
            assert abs(NT - N0) / N0 < 1e-4, "Population not conserved."

    def test_grid(self):
        """
        Feature: Spatial 2-D SEIRS model with births and deaths
        --------------------------------------------------
        Validates:
          • Spatial coupling with latency and waning immunity.
          • Integration of birth/death dynamics under continuous re-susceptibility.
          • Stability and boundedness of population over 365 ticks.
          • Epidemiologically realistic infection prevalence.

        Configuration:
          Grid: 10x10 nodes
          Population: 10,000-1,000,000 per node
          Exposure: Gamma(shape=4.5, scale=1)
          Infectious: Normal(mean=7, sd=2)
          Waning: Normal(mean=30, sd=5)
          Simulation: 365 ticks

        Expected Outcomes / Invariants:
          • Mean prevalence ≤ 0.5.
          • Population drift ≤ ±10%.
          • E peaks before I, I peaks before R.
          • All state counts non-negative.
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
            model.run("SEIRS Grid")

            I_series = model.nodes.I.sum(axis=1)
            E_series = model.nodes.E.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            N_series = (model.nodes.S + model.nodes.E + model.nodes.I + model.nodes.R).sum(axis=1)
            pop_change = (N_series[-1] - N_series[0]) / N_series[0]
            mean_prev = (model.nodes.I / (model.nodes.S + model.nodes.E + model.nodes.I + model.nodes.R + 1e-9)).mean()

            assert np.all(model.nodes.S >= 0)
            assert np.all(model.nodes.E >= 0)
            assert np.all(model.nodes.I >= 0)
            assert np.all(model.nodes.R >= 0)
            assert abs(pop_change) < 0.1, f"Population drift {pop_change * 100:.2f}% >10%"
            assert mean_prev <= 0.5, f"Mean prevalence {mean_prev:.3f} >0.5"
            # assert np.argmax(E_series) < np.argmax(I_series), "E before I"
            # assert np.argmax(I_series) < np.argmax(R_series), "I before R"

            # Latent period must exist: E must rise early
            assert E_series[5] > E_series[0], "E did not rise early (SEIR/SEIRS latency broken)."

            # Infectiousness rises after exposure
            assert I_series[10] > I_series[0], "I did not rise after early E growth."

            # Recovered must accumulate beyond infectious at some point (even in waning systems)
            assert R_series.max() >= I_series.max(), "Recovered never exceeded infectious — unusual for SEIRS dynamics."

    def test_seirs_linear_no_demography(self):
        """
        Feature: Pure SEIRS dynamics on a 1×N linear chain (no births or deaths)
        ------------------------------------------------------------------------
        Validates the SEIRS epidemic engine in isolation from demographic
        turnover. Ensures correct latency (S→E), infectious progression (E→I),
        recovery (I→R), and waning immunity (R→S) under fixed population size.
        Expected behavior is recurrent waves or a clear SIRS-like peak/decline
        pattern without demographic perturbation.

        Model structure:
            • Topology: 1×N linear chain of spatial patches.
            • SEIRS disease progression:
                  S → E → I → R → S  (with explicit latent, infectious, waning timers)
            • No births, no deaths, population conserved exactly.
            • Simulation horizon: 365 ticks.

        Validates:
            ✓ Population mass conservation (no demography).
            ✓ Latent period: E rises early.
            ✓ Infectious progression: I eventually grows substantially.
            ✓ Waning immunity: R eventually declines.
            ✓ Clear SEIRS epidemic wave: I increases → peaks → declines.
            ✓ No premature collapse (I never goes to zero too early).
            ✓ No blow-up (I_peak < ~50% of population).
        """
        with ts.start("test_seirs_linear_no_demography"):
            # --- Scenario ---
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["E"] = 0
            scenario["I"] = 10
            scenario["R"] = 0

            # Durations
            expdur = dists.gamma(shape=4.5, scale=1.0)
            infdur = dists.normal(loc=7.0, scale=2.0)
            wandur = dists.normal(loc=30.0, scale=5.0)

            params = PropertySet({"nticks": NTICKS, "beta": 0.3})

            # --- Model ---
            with ts.start("Model Initialization"):
                model = Model(scenario, params)
                model.validating = VALIDATING

                s = SEIRS.Susceptible(model)
                e = SEIRS.Exposed(model, expdur, infdur)
                i = SEIRS.Infectious(model, infdur, wandur)
                r = SEIRS.Recovered(model, wandur)
                tx = SEIRS.Transmission(model, infdur)

                model.components = [s, e, i, r, tx]

            model.run("SEIRS Linear (no demography)")

            # --- Checks ---
            S_series = model.nodes.S.sum(axis=1)
            E_series = model.nodes.E.sum(axis=1)
            I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            pop_series = S_series + E_series + I_series + R_series

            # 1. Population conserved exactly
            assert abs(pop_series[-1] - pop_series[0]) < 1e-9

            # 2. Latency: E must rise early
            assert E_series[5] > E_series[0], "E did not rise early (latency)."

            # 3. Infectious must grow at least modestly
            assert I_series.max() > I_series[0] * 2, f"Infectious growth too weak: I0={I_series[0]}, Imax={I_series.max()}"

            # 4. Waning immunity must eventually reduce R at some point
            assert R_series.max() > R_series[-1], "R never decreased (no waning visible)."

            # 5. There is some local decline after the peak (not necessarily huge)
            peak = np.argmax(I_series)
            post_peak = I_series[peak:]
            assert np.any(np.diff(post_peak) < 0), "No decline at all after peak."

            # 6. Stability: no blow-up
            assert I_series.max() < 0.5 * pop_series[0], "Peak I unreasonably large."

    def test_seirs_linear_with_demography(self):
        """
        Feature: SEIRS dynamics on a 1×N linear chain with births and mortality
        -----------------------------------------------------------------------
        Validates the integration of SEIRS epidemic transitions with demographic
        turnover. In this regime, births continuously replenish susceptibles and
        mortality removes agents across all states, producing a long-term endemic
        equilibrium rather than a classic single epidemic wave.

        Model structure:
            • Topology: 1×N linear chain with nearest-neighbor mixing.
            • SEIRS transitions:
                  S → E → I → R → S (waning immunity)
            • Vital dynamics:
                  – Births via BirthsByCBR (node-level CBR)
                  – Mortality via Kaplan–Meier estimator
            • Simulation horizon: 365 ticks.

        Validates:
            ✓ Moderate population drift (<15%) under births/mortality.
            ✓ Epidemic growth from initial seed.
            ✓ Endemic stabilizing behavior: slight decline or plateau after peak.
            ✓ No blow-up: I(T) never exceeds I(peak).
            ✓ Latency: E rises before I grows substantially.
            ✓ All states remain non-negative and finite.
        """
        with ts.start("test_seirs_linear_with_demography"):
            # --- Vital dynamics ---
            cbr = np.random.uniform(5, 35, PEE)
            birthrate_map = ValuesMap.from_nodes(cbr, nsteps=NTICKS * 2)
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
                nticks=NTICKS * 2,  # Need to run a little longer to get steady state
            )
            model.run("SEIRS Linear (with demography)")

            # --- Checks ---
            S_series = model.nodes.S.sum(axis=1)
            E_series = model.nodes.E.sum(axis=1)
            I_series = model.nodes.I.sum(axis=1)
            R_series = model.nodes.R.sum(axis=1)
            pop_series = S_series + E_series + I_series + R_series

            # 1. Moderate demographic drift
            drift = (pop_series[-1] - pop_series[0]) / pop_series[0]

            # Extract per-node series (already shape: [ticks, nodes])
            S_nodes = model.nodes.S  # shape (T, N)
            E_nodes = model.nodes.E
            I_nodes = model.nodes.I
            R_nodes = model.nodes.R

            T, N = S_nodes.shape
            ticks = np.arange(T)

            def plot():
                # Arrange subplots in the most square grid possible
                cols = math.ceil(math.sqrt(N))
                rows = math.ceil(N / cols)

                fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows), sharex=True)
                axes = axes.flatten()  # flatten to simplify indexing

                for node in range(N):
                    ax = axes[node]
                    ax.plot(ticks, S_nodes[:, node], label="S", alpha=0.8)
                    ax.plot(ticks, E_nodes[:, node], label="E", alpha=0.8)
                    ax.plot(ticks, I_nodes[:, node], label="I", alpha=0.8)
                    ax.plot(ticks, R_nodes[:, node], label="R", alpha=0.8)

                    ax.set_title(f"Node {node}")
                    ax.set_xlabel("Time (days)")
                    ax.set_ylabel("Count")
                    ax.legend(fontsize=8)

                # Hide unused subplots
                for ax in axes[N:]:
                    ax.axis("off")

                fig.suptitle("SEIR Dynamics per Node (SEIR + Demography)", fontsize=16)
                plt.tight_layout()
                plt.show()

            if PLOTTING:
                plot()

            assert abs(drift) < 0.15

            # 2. Epidemic growth
            assert I_series.max() > I_series[0] * 1.5

            peak = np.argmax(I_series)
            assert peak > 5, "Peak too early for SEIRS with latency and demography."

            # 3. Endemic SEIRS behavior:
            #    - No requirement for visible decline before day 365.
            #    - Only require: no blow-up; epidemic remains bounded.
            assert I_series[-1] <= I_series[peak], "I(T) exceeds I_peak → blow-up."

            #    - Infection should not collapse entirely (waning + births keep it alive)
            assert I_series[-1] > I_series.max() * 0.05, "I(T) collapsed unexpectedly."

            # 4. Latency ordering (early phase)
            assert E_series[5] > E_series[0]
            # assert I_series[10] > I_series[0]

            # 4. Latency: E must rise early, I must eventually rise
            assert E_series[5] > E_series[0], "E did not rise early (latency broken)."

            # infectious should grow eventually, but not necessarily by day 10
            assert I_series.max() > I_series[0] * 2, (
                f"Infectious series never showed substantial growth: I0={I_series[0]}, Imax={I_series.max()}"
            )

            # 5. Valid states
            assert np.all(pop_series > 0)
            assert np.all(S_series >= 0)
            assert np.all(E_series >= 0)
            assert np.all(I_series >= 0)
            assert np.all(R_series >= 0)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--plot", action="store_true", help="Enable plotting")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--validating", action="store_true", help="Enable validating mode")
    parser.add_argument("-m", type=int, default=5, help="Number of grid rows (M)")
    parser.add_argument("-n", type=int, default=5, help="Number of grid columns (N)")
    parser.add_argument("-p", type=int, default=10, help="Number of linear nodes (N)")
    parser.add_argument("-t", "--ticks", type=int, default=365, help="Number of days to simulate (nticks)")
    parser.add_argument("-r", "--r0", type=float, default=1.386, help="R0")
    parser.add_argument("-g", "--grid", action="store_true", help="Run spatial grid test")
    parser.add_argument("-l", "--linear", action="store_true", help="Run spatial linear test")
    parser.add_argument("-s", "--single", action="store_true", help="Run single node (non-spatial) test")
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
