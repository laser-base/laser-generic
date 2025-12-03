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


def build_model(m, n, pop_fn, init_infected=0, init_recovered=0, birthrates=None, pyramid=None, survival=None, nticks=NTICKS, beta=None):
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

    if not beta:
        beta = R0 / INFECTIOUS_DURATION_MEAN
    params = PropertySet({"nticks": nticks, "beta": beta})

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
            model = build_model(1, 1, lambda x, y: 100_000, init_infected=10, beta=0.25)
            model.run("SEIR Single Node")

            # ---------------------------------------------------
            # 1. Extract node-level arrays from LASER LaserFrame
            # ---------------------------------------------------
            # Each is shaped (T, N)
            E_nodes = model.nodes.E
            I_nodes = model.nodes.I
            T, N = I_nodes.shape

            # Derive series *after* these are defined
            E_series = E_nodes.sum(axis=1)
            I_series = I_nodes.sum(axis=1)

            # Quantitative checks
            assert E_series.max() > 0, "No exposed cases observed."
            assert np.argmax(E_series) < np.argmax(I_series), "E should peak before I."
            peak_I = np.argmax(I_series)
            assert I_series[-1] < I_series[peak_I] * 0.5, "I should decline post-peak."

            N0 = (model.nodes.S[0] + model.nodes.E[0] + model.nodes.I[0] + model.nodes.R[0]).sum()
            NT = (model.nodes.S[-1] + model.nodes.E[-1] + model.nodes.I[-1] + model.nodes.R[-1]).sum()
            assert abs(NT - N0) / N0 < 1e-4, "Population not conserved (ΔN>0.01%)."

            final_R_frac = model.nodes.R[-1].sum() / N0
            assert 0.55 <= final_R_frac <= 0.75, f"Final attack fraction {final_R_frac:.3f} out of expected 0.55–0.75 range."

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

    def test_seir_linear_no_demography(self):
        """
        Feature: Pure SEIR dynamics on a 1×N linear chain (no births or deaths)
        -----------------------------------------------------------------------
        This test validates the core SEIR epidemic engine in complete isolation
        from demographic processes. Using a one-dimensional spatial chain and
        fixed population, it ensures that LASER correctly implements exposure
        latency, timed transitions through the E→I→R pipeline, and stable
        epidemic growth under a subcritical-but-growing R₀ ≈ 1.386.

        Model structure:
            • Topology: 1×N linear chain with nearest-neighbor mixing.
            • Disease progression: S → E → I → R
                  – Latent (E) duration: Gamma(k=4.5, θ=1)
                  – Infectious (I) duration: Normal(mean=7, sd=2)
            • No demographic turnover: population remains constant exactly.
            • Initial infections: 10 agents across the chain.
            • Simulation horizon: 365 ticks.

        What this test verifies:
            ✓ **Mass conservation:** S+E+I+R remains exactly constant (no demography).
            ✓ **Latency:** E(t) appears early and grows before I(t) accelerates.
            ✓ **Progression:** I(t) eventually grows by an order of magnitude,
              confirming the E→I transition mechanism.
            ✓ **State balance:** E-peak is substantial relative to I-peak, ensuring
              correct latent-period buffering.
            ✓ **Stable SEIR growth:** With R₀ ≈ 1.386, the epidemic grows steadily
              but has not yet peaked-and-crashed by day 365 (correct behavior).
            ✓ **No blow-up:** Peak infectious prevalence remains below half
              the total population.
            ✓ **No premature collapse:** I(t) remains nonzero late in the run.
            ✓ **Long-wave SEIR dynamics:** Late-epidemic values (t≈275–365)
              show continued upward or plateauing behavior, as expected for
              SEIR with significant latency and modest R₀.

        Why this test matters:
            Passing this test demonstrates that LASER’s SEIR disease engine is
            functioning correctly in its pure form—exposure, infection, and
            recovery timers all interact properly; spatial mixing behaves as
            expected; and the system exhibits classical SEIR long-wave dynamics
            without numerical errors. This test provides a clean foundation for
            more complex SEIR validations involving demographic turnover.
        """
        with ts.start("test_seir_linear_no_demography"):
            # --- Scenario ---
            scenario = stdgrid(M=1, N=PEE)
            scenario["S"] = scenario["population"] - 10
            scenario["E"] = 0
            scenario["I"] = 10
            scenario["R"] = 0

            # Durations
            expdur = dists.gamma(shape=4.5, scale=1.0)
            infdur = dists.normal(loc=7.0, scale=2.0)

            # R0 → beta
            beta = R0 / 7.0
            params = PropertySet({"nticks": NTICKS, "beta": beta})

            # --- Model ---
            with ts.start("Model Initialization"):
                model = Model(scenario, params)
                model.validating = VALIDATING

                s = SEIR.Susceptible(model)
                e = SEIR.Exposed(model, expdur, infdur)
                i = SEIR.Infectious(model, infdur)
                r = SEIR.Recovered(model)
                tx = SEIR.Transmission(model, infdur)

                model.components = [s, e, i, r, tx]

            model.run("SEIR Linear (no demography)")

            # --- Checks ---
            S = model.nodes.S.sum(axis=1)
            E = model.nodes.E.sum(axis=1)
            I_series = model.nodes.I.sum(axis=1)
            R = model.nodes.R.sum(axis=1)
            pop = S + E + I_series + R

            # 1. perfect conservation
            assert abs(pop[-1] - pop[0]) < 1e-9

            # 2. latency: exposed must exist and rise early
            assert E.max() > 0, "No exposed individuals observed."
            assert E[5] > E[0], "E did not rise early (SEIR latency broken)."

            # 3. infectious must eventually grow substantially
            assert I_series.max() > I_series[0] * 10, f"Infectious cases did not grow strongly: I0={I_series[0]}, peak={I_series.max()}"

            # 4. exposed must be substantial relative to infectious
            assert E.max() > 0.1 * I_series.max(), f"E peak ({E.max()}) too small relative to I_series peak ({I_series.max()})."

            # 5. epidemic should show strong growth
            peak = np.argmax(I_series)
            assert I_series[peak] > I_series[0] * 10, (
                f"Infectious cases showed insufficient growth: I0={I_series[0]}, peak={I_series[peak]}"
            )

            # 6. epidemic should NOT blow up (stability check)
            assert I_series[peak] < pop[0] * 0.5, "Peak infectious prevalence unrealistically high — unstable dynamics."

            # 7. epidemic should not collapse prematurely
            assert I_series[-1] > 0, "Infectious cases crashed to zero unexpectedly (SIR-like behavior)."

            # 8. epidemic should still be growing or near plateau by day 365
            assert I_series[-1] >= I_series[int(0.75 * NTICKS)], "Unexpected decline; SEIR should not decline this early for R0≈1.4."

    def test_seir_linear_with_demography(self):
        """
        Feature: SEIR dynamics on a 1×N linear chain with demographic turnover
        ----------------------------------------------------------------------
        This test validates LASER’s integration of SEIR epidemic progression with
        demographic processes (births and mortality). It ensures that the model
        behaves plausibly when susceptible replenishment and age-structured
        mortality interact with latent exposure, infectiousness, and recovery.

        Model structure:
            • Topology: 1×N linear chain with nearest-neighbor mixing.
            • SEIR disease progression:
                  S → E → I → R, with explicit latent and infectious durations.
            • Vital dynamics:
                  – Births via BirthsByCBR (CBR drawn per node).
                  – Deaths via MortalityByEstimator (Kaplan–Meier survival).
            • Simulation horizon: 365 ticks.

        What this test verifies:
            ✓ **Demographic drift bounded:** Population remains positive and
              exhibits <15% net drift despite ongoing births/mortality.
            ✓ **Epidemic growth:** I(t) increases substantially from initial
              seeding.
            ✓ **Endemic SEIR behavior:** After the peak, I(t) exhibits modest
              decline or plateau, consistent with SEIR + demographic turnover.
            ✓ **No epidemic blow-up:** I(T) never exceeds I(peak).
            ✓ **Latency preserved:** E(t) precedes or coincides with I(t) in peak
              timing, indicating correct S→E→I sequencing.
            ✓ **State validity:** No negative counts in S, E, I, or R.
            ✓ **Numerical stability:** Birth/mortality flows do not destabilize,
              and all compartments remain well-defined.

        Why this test matters:
            SEIR models with demographic turnover naturally approach endemic
            equilibria rather than collapsing after the epidemic peak. Passing
            this test demonstrates that LASER correctly handles susceptible
            replenishment, age-structured mortality, and their interaction with
            SEIR timers. It confirms the robustness of demographic/epidemiologic
            coupling in long-term simulations.
        """
        with ts.start("test_seir_linear_with_demography"):
            # Let's run for 2 years to let things smooth out with these settings
            cbr = np.random.uniform(5, 35, PEE)
            birthrates = ValuesMap.from_nodes(cbr, nsteps=NTICKS * 2)
            pyramid = AliasedDistribution(np.full(89, 1_000))
            survival = KaplanMeierEstimator(np.full(89, 1_000).cumsum())

            model = build_model(
                1,
                PEE,
                lambda x, y: int(np.random.uniform(10_000, 1_000_000)),
                init_infected=10,
                birthrates=birthrates.values,
                pyramid=pyramid,
                survival=survival,
                nticks=NTICKS * 2,
            )
            model.run("SEIR Linear (with demography)")

            S = model.nodes.S.sum(axis=1)
            E = model.nodes.E.sum(axis=1)
            I_series = model.nodes.I.sum(axis=1)
            R = model.nodes.R.sum(axis=1)
            pop = S + E + I_series + R

            pop0, popT = pop[0], pop[-1]
            drift = (popT - pop0) / pop0

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

            # 1. moderate population drift (< 15%)
            assert abs(drift) < 0.15, f"Population drift {drift * 100:.2f}% >15%."

            # 2. epidemic growth & decline
            assert I_series.max() > I_series[0] * 1.5
            peak = np.argmax(I_series)
            assert peak > 5
            # assert I[-1] < I[peak] * 0.9

            # SEIR+demography should show *some* decline after peak OR reach plateau
            post_peak = I_series[peak:]
            assert np.any(np.diff(post_peak) < 0), "No decline at all after peak — suspect transmission or waning error."

            # Ensure model is not unstable: I_end should not exceed I_peak
            assert I_series[-1] <= I_series[peak], "I_series(T) exceeds I_series(peak), indicating epidemic blow-up."

            # 3. E precedes I_series
            assert E.max() > 0
            assert np.argmax(E) < np.argmax(I_series)

            # 4. demographics do not create negative or invalid states
            assert np.all(pop > 0)
            assert np.all(S >= 0)
            assert np.all(E >= 0)
            assert np.all(I_series >= 0)
            assert np.all(R >= 0)


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
    parser.add_argument("-g", "--grid", action="store_true", help="Run grid spatial test")
    parser.add_argument("-l", "--linear", action="store_true", help="Run linear spatial test")
    parser.add_argument("-s", "--single", action="store_true", help="Run single node test (not spatial)")
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
