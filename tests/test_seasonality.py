"""
Comprehensive tests for seasonal forcing in transmission components.

Tests three transmission types (SI, SIR, SEIR) with three seasonality
conditions (none, attenuated, amplified) for a total of 9 test scenarios.
"""

import unittest

import laser.core.distributions as dists
import numpy as np
from laser.core import PropertySet
from laser.core.random import seed as set_seed
from scipy.special import lambertw

from laser.generic import SEIR
from laser.generic import SI
from laser.generic import SIR
from laser.generic import Model
from laser.generic.components import TransmissionSE
from laser.generic.components import TransmissionSI
from laser.generic.components import TransmissionSIx
from laser.generic.utils import ValuesMap

try:
    from tests.utils import stdgrid
except ImportError:
    from utils import stdgrid

# Shared test parameters
NTICKS = 730  # 2 years to observe seasonal patterns
SEED = 271828
POPULATION = 100_000


def create_seasonality_valuesmap(multiplier: float, nnodes: int, nticks: int) -> ValuesMap:
    """Create constant seasonality forcing."""
    return ValuesMap.from_timeseries(np.full(nticks, multiplier, dtype=np.float32), nnodes)


def calculate_infection_rate(model, start=1, end=40):
    """
    Calculate average infection rate from SI model during growth phase.
    Rate = dI / (S * dt)
    Uses early time window (1-40 days) to capture exponential growth.
    """
    infected_series = model.nodes.I[start:end].sum(axis=1)
    susceptible_series = model.nodes.S[start:end].sum(axis=1)
    delta_infected = np.diff(infected_series)
    S_avg = (susceptible_series[:-1] + susceptible_series[1:]) / 2

    # Avoid division by zero
    valid_mask = S_avg > 0
    if not np.any(valid_mask):
        return 0.0

    rates = delta_infected[valid_mask] / S_avg[valid_mask]
    return np.mean(rates)


def calculate_attack_fraction(beta, inf_mean, pop, init_inf):
    """Calculate final attack rate using Kermack-McKendrick formula."""
    R0 = beta * inf_mean
    S0 = (pop - init_inf) / pop
    S_inf = -1 / R0 * lambertw(-R0 * S0 * np.exp(-R0)).real
    return float(1 - S_inf)


class TestSeasonalForcing(unittest.TestCase):
    """
    Comprehensive tests for seasonal forcing in transmission components.

    Tests three transmission types (SI, SIR, SEIR) with three seasonality
    conditions (none, attenuated, amplified) for a total of 9 test scenarios.
    """

    def setUp(self):
        """Set random seed for reproducibility."""
        set_seed(SEED)
        return

    # ========================================================================
    # SI Model Tests (TransmissionSIx)
    # ========================================================================

    def test_si_no_seasonality(self):
        """
        Validate baseline SI model transmission dynamics without seasonal forcing.

        WHAT IS VALIDATED:
        - Epidemic growth occurs (final infections > initial infections)
        - Baseline infection rate is positive during exponential growth phase

        SCENARIO SETUP:
        - Single node with 100,000 population
        - 1,000 initial infections, 99,000 susceptible
        - Beta = 0.5, no seasonality multiplier (baseline transmission)
        - 730-day simulation (2 years)

        FAILURE MEANING:
        If this test fails, the baseline SI transmission model is broken. The epidemic
        should always grow in a fully susceptible population with positive beta. Failure
        indicates a fundamental problem with TransmissionSIx component or SI model logic.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 1000  # Start with 1000 infected
        scenario["I"] = 1000

        params = PropertySet({"nticks": NTICKS, "beta": 0.5})
        model = Model(scenario, params)

        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=None),
        ]

        model.run("SI No Seasonality")

        # Epidemic should grow
        initial_I = model.nodes.I[0].sum()
        final_I = model.nodes.I[-1].sum()
        assert final_I > initial_I, "SI epidemic should grow over time"

        # Calculate baseline infection rate for comparison
        baseline_rate = calculate_infection_rate(model)
        assert baseline_rate > 0, "Baseline infection rate should be positive"

        return

    def test_si_attenuated_seasonality(self):
        """
        Validate that 0.5× seasonal forcing correctly reduces SI transmission rates.

        WHAT IS VALIDATED:
        - Infection rate is approximately 50% of baseline (within 30-70% range)
        - Time to 50% prevalence is delayed compared to baseline
        - Force of infection (FOI) is approximately 50% of baseline during growth phase

        SCENARIO SETUP:
        - Single node with 100,000 population, 1,000 initial infections
        - Beta = 0.5, seasonality = 0.5× (constant attenuation)
        - Compares against baseline run without seasonality
        - Measures infection rate during days 1-40 (exponential growth window)

        FAILURE MEANING:
        If this test fails, seasonal forcing is not properly modulating transmission. The 0.5×
        multiplier should halve the effective transmission rate, slowing the epidemic proportionally.
        Failure suggests the seasonality parameter is not being applied to the force of infection
        calculation in TransmissionSIx, or the multiplier is being applied incorrectly.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 1000
        scenario["I"] = 1000

        params = PropertySet({"nticks": NTICKS, "beta": 0.5})
        model = Model(scenario, params)

        seasonality = create_seasonality_valuesmap(0.5, nnodes=1, nticks=NTICKS)

        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=seasonality),
        ]

        model.run("SI Attenuated Seasonality")

        # Get baseline for comparison
        baseline_model = self._create_si_baseline()
        baseline_rate = calculate_infection_rate(baseline_model)

        # Calculate attenuated metrics
        attenuated_rate = calculate_infection_rate(model)

        # Infection rate should be approximately 50% of baseline (+/- 20% tolerance)
        rate_ratio = attenuated_rate / baseline_rate
        assert 0.3 < rate_ratio < 0.7, f"Attenuated infection rate ratio {rate_ratio:.3f} should be ~0.5 (+/- 0.2)"

        # Epidemic should reach 50% prevalence later than baseline
        baseline_half_time = np.argmax(baseline_model.nodes.I[:, 0] >= POPULATION / 2)
        attenuated_half_time = np.argmax(model.nodes.I[:, 0] >= POPULATION / 2)
        assert attenuated_half_time > baseline_half_time, (
            f"Attenuated epidemic should reach 50% prevalence later (day {attenuated_half_time}) than baseline (day {baseline_half_time})"
        )

        # FOI should be consistently lower during growth phase
        baseline_foi = baseline_model.nodes.forces[:, 0]
        attenuated_foi = model.nodes.forces[:, 0]
        mean_foi_ratio = np.mean(attenuated_foi[1:40]) / np.mean(baseline_foi[1:40])
        assert 0.3 < mean_foi_ratio < 0.7, f"Mean FOI ratio {mean_foi_ratio:.3f} should be ~0.5"

        return

    def test_si_amplified_seasonality(self):
        """
        Validate that 2.0× seasonal forcing correctly amplifies SI transmission rates.

        WHAT IS VALIDATED:
        - Infection rate is approximately 200% of baseline (within 140-260% range)
        - Time to 50% prevalence is accelerated compared to baseline
        - Force of infection (FOI) is approximately 200% of baseline during growth phase

        SCENARIO SETUP:
        - Single node with 100,000 population, 1,000 initial infections
        - Beta = 0.5, seasonality = 2.0× (constant amplification)
        - Compares against baseline run without seasonality
        - Measures infection rate during days 1-40 (exponential growth window)

        FAILURE MEANING:
        If this test fails, seasonal forcing amplification is not working correctly. The 2.0×
        multiplier should double the effective transmission rate, accelerating the epidemic
        proportionally. Failure indicates the seasonality parameter is not properly scaling
        up transmission in TransmissionSIx, possibly due to capping or incorrect multiplication.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 1000
        scenario["I"] = 1000

        params = PropertySet({"nticks": NTICKS, "beta": 0.5})
        model = Model(scenario, params)

        seasonality = create_seasonality_valuesmap(2.0, nnodes=1, nticks=NTICKS)

        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=seasonality),
        ]

        model.run("SI Amplified Seasonality")

        # Get baseline for comparison
        baseline_model = self._create_si_baseline()
        baseline_rate = calculate_infection_rate(baseline_model)

        # Calculate amplified metrics
        amplified_rate = calculate_infection_rate(model)

        # Infection rate should be approximately 200% of baseline (+/- 30% tolerance)
        rate_ratio = amplified_rate / baseline_rate
        assert 1.4 < rate_ratio < 2.6, f"Amplified infection rate ratio {rate_ratio:.3f} should be ~2.0 (+/- 0.6)"

        # Epidemic should reach 50% prevalence earlier than baseline
        baseline_half_time = np.argmax(baseline_model.nodes.I[:, 0] >= POPULATION / 2)
        amplified_half_time = np.argmax(model.nodes.I[:, 0] >= POPULATION / 2)
        assert amplified_half_time < baseline_half_time, (
            f"Amplified epidemic should reach 50% prevalence earlier (day {amplified_half_time}) than baseline (day {baseline_half_time})"
        )

        # FOI should be consistently higher during growth phase
        baseline_foi = baseline_model.nodes.forces[:, 0]
        amplified_foi = model.nodes.forces[:, 0]
        mean_foi_ratio = np.mean(amplified_foi[1:40]) / np.mean(baseline_foi[1:40])
        assert 1.4 < mean_foi_ratio < 2.6, f"Mean FOI ratio {mean_foi_ratio:.3f} should be ~2.0"

        return

    # ========================================================================
    # SIR Model Tests (TransmissionSI)
    # ========================================================================

    def test_sir_no_seasonality(self):
        """
        Validate baseline SIR model with Kermack-McKendrick theoretical prediction.

        WHAT IS VALIDATED:
        - Epidemic curve shows characteristic rise and fall with peak in simulation window
        - Final attack rate matches K-M analytical prediction (within 10% tolerance)
        - Population is conserved (S + I + R = constant)

        SCENARIO SETUP:
        - Single node with 100,000 population, 500 initial infections
        - Beta = 0.25, infectious period mean = 7 days, giving R₀ ≈ 1.75
        - No seasonality (baseline transmission)
        - K-M formula predicts final size based on R₀ and initial conditions

        FAILURE MEANING:
        If this test fails, the SIR model implementation deviates from epidemic theory. The
        Kermack-McKendrick equation provides an analytical solution for final attack rate
        given R₀. Failure suggests errors in recovery timing (infdurdist), transmission
        calculation (TransmissionSI), or state transitions between S/I/R compartments.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 500  # Start with 500 infected
        scenario["I"] = 500
        scenario["R"] = 0

        beta = 0.25
        inf_mean = 7.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)

        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            SIR.Susceptible(model),
            SIR.Infectious(model, infdurdist),
            SIR.Recovered(model),
            TransmissionSI(model, infdurdist, seasonality=None),
        ]

        model.run("SIR No Seasonality")

        # Epidemic curve should show rise and fall
        I_series = model.nodes.I[:, 0]
        peak_time = np.argmax(I_series)
        assert 10 < peak_time < NTICKS - 50, "Peak should occur during simulation with time to resolve"

        # Final attack rate should match Kermack-McKendrick (+/- 10% tolerance)
        expected_af = calculate_attack_fraction(beta, inf_mean, POPULATION, 500)
        actual_af = model.nodes.R[-1].sum() / POPULATION
        percent_diff = abs((actual_af - expected_af) / expected_af * 100)

        assert percent_diff < 10.0, f"Attack fraction {actual_af:.4f} deviates by {percent_diff:.1f}% from K-M prediction {expected_af:.4f}"

        # Population should be conserved
        total_pop = model.nodes.S + model.nodes.I + model.nodes.R
        assert np.all(total_pop == POPULATION), "Population should be conserved"

        return

    def test_sir_attenuated_seasonality(self):
        """
        Validate SIR epidemic die-out when seasonal forcing reduces R₀ below 1.

        WHAT IS VALIDATED:
        - Final attack rate is minimal (<15%) indicating epidemic extinction
        - Peak infections occur early (<100 days) and remain low (<3,500)
        - Epidemic dies out rather than spreading through population

        SCENARIO SETUP:
        - Single node with 100,000 population, 2,000 initial infections
        - Beta = 0.25, infectious period = 7 days, seasonality = 0.5×
        - Effective R₀ = 0.5 × 0.25 × 7 ≈ 0.875 < 1 (subcritical)
        - Larger initial infected to observe dynamics before die-out

        FAILURE MEANING:
        If this test fails, seasonality is not reducing transmission below the epidemic threshold.
        With R₀ < 1, each infected person infects <1 other on average, so the epidemic must die out.
        Failure indicates seasonality isn't properly attenuating transmission, or recovery rates
        are incorrect. An attack rate >15% suggests R₀ is still above 1 despite the multiplier.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        # Use larger initial infected to slow die-out
        scenario["S"] = scenario.population - 2000
        scenario["I"] = 2000
        scenario["R"] = 0

        beta = 0.25
        inf_mean = 7.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)

        # With 0.5x seasonality: effective beta = 0.125, R0 approx 0.875 < 1
        seasonality = create_seasonality_valuesmap(0.5, nnodes=1, nticks=NTICKS)
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            SIR.Susceptible(model),
            SIR.Infectious(model, infdurdist),
            SIR.Recovered(model),
            TransmissionSI(model, infdurdist, seasonality=seasonality),
        ]

        model.run("SIR Attenuated Seasonality")

        # Since R0 < 1, epidemic should die out
        final_attack_rate = model.nodes.R[-1].sum() / POPULATION

        # Attack rate should be minimal (much less than baseline ~60%)
        assert final_attack_rate < 0.15, f"With R0 < 1, attack rate {final_attack_rate:.4f} should be < 0.15"

        # I compartment should decrease monotonically after initial period
        I_series = model.nodes.I[:, 0]
        peak_I = np.max(I_series)
        peak_time = np.argmax(I_series)

        # Peak should be early (within first 100 days) and relatively low
        assert peak_time < 100, f"With R0 < 1, peak should occur early (at day {peak_time})"
        assert peak_I < 3500, f"With R0 < 1, peak infections {peak_I} should be low"

        return

    def test_sir_amplified_seasonality(self):
        """
        Validate SIR large outbreak dynamics with amplified seasonal forcing (R₀ ≈ 3.5).

        WHAT IS VALIDATED:
        - Final attack rate matches K-M prediction for R₀ ≈ 3.5 (within 10%)
        - Attack rate exceeds baseline (no seasonality) scenario
        - Peak infections are substantial (>10,000) and occur early (<200 days)

        SCENARIO SETUP:
        - Single node with 100,000 population, 500 initial infections
        - Beta = 0.25, infectious period = 7 days, seasonality = 2.0×
        - Effective R₀ = 2.0 × 0.25 × 7 ≈ 3.5 (supercritical)
        - K-M formula predicts ~95% final attack rate for R₀ = 3.5

        FAILURE MEANING:
        If this test fails, seasonal amplification is not properly scaling transmission upward.
        With R₀ = 3.5, we expect a severe, fast outbreak affecting most of the population. Failure
        suggests the 2.0× multiplier isn't being applied correctly, peak timing/magnitude indicates
        slower-than-expected transmission, or K-M deviation suggests incorrect compartment dynamics.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 500
        scenario["I"] = 500
        scenario["R"] = 0

        beta = 0.25
        inf_mean = 7.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)

        # With 2.0x seasonality: effective beta = 0.5, R0 approx 3.5
        seasonality = create_seasonality_valuesmap(2.0, nnodes=1, nticks=NTICKS)
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            SIR.Susceptible(model),
            SIR.Infectious(model, infdurdist),
            SIR.Recovered(model),
            TransmissionSI(model, infdurdist, seasonality=seasonality),
        ]

        model.run("SIR Amplified Seasonality")

        # Calculate expected attack fraction with effective R0
        effective_beta = beta * 2.0
        expected_af = calculate_attack_fraction(effective_beta, inf_mean, POPULATION, 500)
        actual_af = model.nodes.R[-1].sum() / POPULATION

        # Should match K-M prediction (+/- 10% tolerance)
        percent_diff = abs((actual_af - expected_af) / expected_af * 100)
        assert percent_diff < 10.0, (
            f"Attack fraction {actual_af:.4f} deviates by {percent_diff:.1f}% from K-M prediction {expected_af:.4f} for R0 approx 3.5"
        )

        # Should have larger outbreak than baseline
        baseline_af = calculate_attack_fraction(beta, inf_mean, POPULATION, 500)
        assert actual_af > baseline_af, f"Amplified attack rate {actual_af:.4f} should exceed baseline {baseline_af:.4f}"

        # Peak should be higher and earlier
        I_series = model.nodes.I[:, 0]
        peak_I = np.max(I_series)
        peak_time = np.argmax(I_series)

        assert peak_I > 10000, f"With R0 approx 3.5, peak {peak_I} should be substantial"
        assert peak_time < 200, f"With R0 approx 3.5, peak should occur early (at day {peak_time})"

        return

    # ========================================================================
    # SEIR Model Tests (TransmissionSE)
    # ========================================================================

    def test_seir_no_seasonality(self):
        """
        Validate baseline SEIR model with exposed compartment dynamics and K-M comparison.

        WHAT IS VALIDATED:
        - Exposed (E) compartment peaks before Infectious (I) compartment
        - Final attack rate approximately matches K-M prediction (within 20% tolerance)
        - Population is conserved (S + E + I + R = constant)

        SCENARIO SETUP:
        - Single node with 100,000 population, 500 initial infections
        - Beta = 0.25, exposed period = 5 days, infectious period = 7 days
        - R₀ ≈ 1.75 (same as SIR for comparison)
        - K-M prediction is approximate since K-M theory assumes SIR, not SEIR

        FAILURE MEANING:
        If this test fails, the SEIR exposed compartment dynamics are incorrect. E must peak
        before I because people progress E→I→R. Failure indicates incorrect transition timing
        (expdurdist/infdurdist), wrong compartment updates, or transmission applied to wrong
        state. K-M tolerance is 20% (vs 10% for SIR) since K-M doesn't account for E compartment.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 500
        scenario["E"] = 0
        scenario["I"] = 500
        scenario["R"] = 0

        beta = 0.25
        exp_mean = 5.0
        inf_mean = 7.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)

        expdurdist = dists.normal(loc=exp_mean, scale=1.0)
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            SEIR.Susceptible(model),
            SEIR.Exposed(model, expdurdist, infdurdist),
            SEIR.Infectious(model, infdurdist),
            SEIR.Recovered(model),
            TransmissionSE(model, expdurdist, seasonality=None),
        ]

        model.run("SEIR No Seasonality")

        # E should peak before I
        E_series = model.nodes.E[:, 0]
        I_series = model.nodes.I[:, 0]
        E_peak_time = np.argmax(E_series)
        I_peak_time = np.argmax(I_series)

        assert E_peak_time < I_peak_time, f"E should peak (day {E_peak_time}) before I (day {I_peak_time})"

        # Final attack rate should be reasonable (K-M is approximate for SEIR due to exposed state)
        # Allow 20% tolerance since K-M theory doesn't account for exposed compartment
        expected_af = calculate_attack_fraction(beta, inf_mean, POPULATION, 500)
        actual_af = model.nodes.R[-1].sum() / POPULATION
        percent_diff = abs((actual_af - expected_af) / expected_af * 100)

        assert percent_diff < 20.0, f"Attack fraction {actual_af:.4f} deviates by {percent_diff:.1f}% from K-M prediction {expected_af:.4f}"

        # Population should be conserved
        total_pop = model.nodes.S + model.nodes.E + model.nodes.I + model.nodes.R
        assert np.all(total_pop == POPULATION), "Population should be conserved"

        return

    def test_seir_attenuated_seasonality(self):
        """
        Validate SEIR epidemic extinction when seasonal forcing reduces R₀ below 1.

        WHAT IS VALIDATED:
        - E and I compartments approach zero (<50 each) at simulation end
        - Final attack rate is minimal (<15%) indicating die-out
        - Peak E and I are low (<3,000) indicating failed outbreak

        SCENARIO SETUP:
        - Single node with 100,000 population, 2,000 initial infections
        - Beta = 0.25, exposed = 5d, infectious = 7d, seasonality = 0.5×
        - Effective R₀ ≈ 0.875 < 1 (subcritical)
        - Larger initial infections to observe dynamics before extinction

        FAILURE MEANING:
        If this test fails, the attenuated seasonality isn't reducing transmission below threshold.
        With R₀ < 1, the epidemic cannot sustain itself and must die out. Persistence of E or I
        compartments at simulation end indicates ongoing transmission despite R₀ < 1. High attack
        rates suggest effective R₀ is still >1, meaning seasonality isn't properly applied in
        TransmissionSE or the multiplier value is not respected.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        # Larger initial infected to observe dynamics before die-out
        scenario["S"] = scenario.population - 2000
        scenario["E"] = 0
        scenario["I"] = 2000
        scenario["R"] = 0

        beta = 0.25
        exp_mean = 5.0
        inf_mean = 7.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)

        # With 0.5x seasonality: effective R0 approx 0.875 < 1
        seasonality = create_seasonality_valuesmap(0.5, nnodes=1, nticks=NTICKS)
        expdurdist = dists.normal(loc=exp_mean, scale=1.0)
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            SEIR.Susceptible(model),
            SEIR.Exposed(model, expdurdist, infdurdist),
            SEIR.Infectious(model, infdurdist),
            SEIR.Recovered(model),
            TransmissionSE(model, expdurdist, seasonality=seasonality),
        ]

        model.run("SEIR Attenuated Seasonality")

        # Epidemic should die out
        final_E = model.nodes.E[-1].sum()
        final_I = model.nodes.I[-1].sum()

        # E and I should be near zero at end
        assert final_E < 50, f"E compartment should be near zero (got {final_E})"
        assert final_I < 50, f"I compartment should be near zero (got {final_I})"

        # Minimal final attack rate
        final_attack_rate = model.nodes.R[-1].sum() / POPULATION
        assert final_attack_rate < 0.15, f"With R0 < 1, attack rate {final_attack_rate:.4f} should be < 0.15"

        # Peak prevalence should be low
        E_series = model.nodes.E[:, 0]
        I_series = model.nodes.I[:, 0]
        peak_E = np.max(E_series)
        peak_I = np.max(I_series)

        assert peak_E < 3000, f"With R0 < 1, peak E {peak_E} should be low"
        assert peak_I < 3000, f"With R0 < 1, peak I {peak_I} should be low"

        return

    def test_seir_amplified_seasonality(self):
        """
        Validate SEIR large outbreak with amplified seasonal forcing (R₀ ≈ 3.5).

        WHAT IS VALIDATED:
        - E peaks before I even with amplified transmission
        - Final attack rate matches K-M prediction for R₀ ≈ 3.5 (within 20%)
        - Attack rate exceeds baseline scenario
        - Peak E and I are substantial (>5,000 each)

        SCENARIO SETUP:
        - Single node with 100,000 population, 500 initial infections
        - Beta = 0.25, exposed = 5d, infectious = 7d, seasonality = 2.0×
        - Effective R₀ = 2.0 × 0.25 × 7 ≈ 3.5 (severe outbreak)
        - K-M comparison has 20% tolerance due to E compartment

        FAILURE MEANING:
        If this test fails, seasonal amplification isn't properly scaling SEIR transmission. The 2.0×
        multiplier should create a severe, rapid outbreak with high attack rates. E-before-I failure
        indicates broken compartment progression. Low peak values suggest transmission isn't being
        amplified as expected. K-M deviation suggests the effective R₀ achieved doesn't match the
        theoretical R₀ = seasonality × beta × infectious_period.
        """

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 500
        scenario["E"] = 0
        scenario["I"] = 500
        scenario["R"] = 0

        beta = 0.25
        exp_mean = 5.0
        inf_mean = 7.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)

        # With 2.0x seasonality: effective R0 approx 3.5
        seasonality = create_seasonality_valuesmap(2.0, nnodes=1, nticks=NTICKS)
        expdurdist = dists.normal(loc=exp_mean, scale=1.0)
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            SEIR.Susceptible(model),
            SEIR.Exposed(model, expdurdist, infdurdist),
            SEIR.Infectious(model, infdurdist),
            SEIR.Recovered(model),
            TransmissionSE(model, expdurdist, seasonality=seasonality),
        ]

        model.run("SEIR Amplified Seasonality")

        # E should still peak before I
        E_series = model.nodes.E[:, 0]
        I_series = model.nodes.I[:, 0]
        E_peak_time = np.argmax(E_series)
        I_peak_time = np.argmax(I_series)

        assert E_peak_time < I_peak_time, "E should peak before I even with amplification"

        # Calculate expected attack fraction with effective R0
        effective_beta = beta * 2.0
        expected_af = calculate_attack_fraction(effective_beta, inf_mean, POPULATION, 500)
        actual_af = model.nodes.R[-1].sum() / POPULATION

        # Should be close to K-M prediction (allow 20% tolerance for SEIR)
        percent_diff = abs((actual_af - expected_af) / expected_af * 100)
        assert percent_diff < 20.0, (
            f"Attack fraction {actual_af:.4f} deviates by {percent_diff:.1f}% from K-M prediction {expected_af:.4f} for R0 approx 3.5"
        )

        # Should have larger outbreak than baseline
        baseline_af = calculate_attack_fraction(beta, inf_mean, POPULATION, 500)
        assert actual_af > baseline_af, f"Amplified attack rate {actual_af:.4f} should exceed baseline {baseline_af:.4f}"

        # Peak prevalence should be substantial
        peak_E = np.max(E_series)
        peak_I = np.max(I_series)

        assert peak_E > 5000, f"With R0 approx 3.5, peak E {peak_E} should be substantial"
        assert peak_I > 5000, f"With R0 approx 3.5, peak I {peak_I} should be substantial"

        return

    # ========================================================================
    # Spatially Varying Seasonality Tests
    # ========================================================================

    def test_si_spatial_seasonality(self):
        """
        Validate spatially independent seasonal forcing across disconnected nodes (SI model).

        WHAT IS VALIDATED:
        - Time to 50% prevalence follows ordering: amplified < baseline < attenuated
        - Three independent nodes with different seasonality respond independently
        - No spatial coupling affects results (gravity_k=0 disables spatial transmission)

        SCENARIO SETUP:
        - Three nodes (1×3 grid), each with 100,000 population, 1,000 initial infections
        - Beta = 0.5, seasonality = [0.5×, 1.0×, 2.0×] for nodes [0, 1, 2]
        - gravity_k = 0.0 (spatial transmission disabled - nodes are isolated)
        - ValuesMap.from_nodes() assigns different constant seasonality per node

        FAILURE MEANING:
        If this test fails, spatial seasonality isn't working correctly. Each node should behave
        like an independent SI model with its own seasonal multiplier. Failure indicates seasonality
        isn't being applied per-node (maybe using global value), spatial transmission isn't actually
        disabled (nodes affecting each other), or ValuesMap.from_nodes() isn't correctly mapping
        values. Incorrect timing order suggests wrong seasonality applied to wrong node.
        """

        # Create 3-node scenario with no spatial transmission
        scenario = stdgrid(M=1, N=3, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 1000
        scenario["I"] = 1000

        # Disable spatial transmission by setting gravity_k to 0
        params = PropertySet({"nticks": NTICKS, "beta": 0.5, "gravity_k": 0.0})
        model = Model(scenario, params)

        # Create spatially varying seasonality: [0.5, 1.0, 2.0]
        seasonality_values = np.array([0.5, 1.0, 2.0], dtype=np.float32)
        seasonality = ValuesMap.from_nodes(seasonality_values, NTICKS)

        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=seasonality),
        ]

        model.run("SI Spatial Seasonality")

        # Calculate time to reach 50% prevalence for each node
        half_pop = POPULATION / 2
        half_times = []
        for node_idx in range(3):
            I_series = model.nodes.I[:, node_idx]
            half_time = np.argmax(I_series >= half_pop)
            half_times.append(half_time)

        attenuated_time, baseline_time, amplified_time = half_times

        # Attenuated node should reach 50% later than baseline
        assert attenuated_time > baseline_time, (
            f"Attenuated node (0.5x) should reach 50% later (day {attenuated_time}) than baseline (day {baseline_time})"
        )

        # Amplified node should reach 50% earlier than baseline
        assert amplified_time < baseline_time, (
            f"Amplified node (2.0x) should reach 50% earlier (day {amplified_time}) than baseline (day {baseline_time})"
        )

        # Verify ordering: amplified < baseline < attenuated
        assert amplified_time < baseline_time < attenuated_time, (
            f"Time ordering should be amplified ({amplified_time}) < baseline ({baseline_time}) < attenuated ({attenuated_time})"
        )

        return

    def test_sir_spatial_seasonality(self):
        """
        Validate spatially independent seasonal forcing across disconnected nodes (SIR model).

        WHAT IS VALIDATED:
        - Amplified node (R₀=3.5) peaks before baseline node (R₀=1.75)
        - Attack rates follow ordering: amplified > baseline > attenuated
        - Attenuated node (R₀<1) has minimal attack rate from die-out

        SCENARIO SETUP:
        - Three nodes (1×3 grid), each with 100,000 population, 500 initial infections
        - Beta = 0.25, inf_period = 7d, seasonality = [0.5×, 1.0×, 2.0×]
        - Effective R₀ = [0.875, 1.75, 3.5] - subcritical, moderate, severe
        - gravity_k = 0.0 (nodes isolated)

        FAILURE MEANING:
        If this test fails, spatial seasonality isn't correctly creating different epidemic outcomes
        in different locations. Node 0 should see die-out (R₀<1), node 1 moderate outbreak (R₀≈1.75),
        node 2 severe outbreak (R₀≈3.5). Failure indicates per-node seasonality not applied, spatial
        transmission still coupling nodes despite gravity_k=0, or ValuesMap.from_nodes() not working.
        Wrong attack rate ordering suggests seasonality values assigned to wrong nodes.
        """

        # Create 3-node scenario with no spatial transmission
        scenario = stdgrid(M=1, N=3, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 500
        scenario["I"] = 500
        scenario["R"] = 0

        beta = 0.25
        inf_mean = 7.0

        # Disable spatial transmission by setting gravity_k to 0
        params = PropertySet({"nticks": NTICKS, "beta": beta, "gravity_k": 0.0})
        model = Model(scenario, params)

        # Create spatially varying seasonality: [0.5, 1.0, 2.0]
        seasonality_values = np.array([0.5, 1.0, 2.0], dtype=np.float32)
        seasonality = ValuesMap.from_nodes(seasonality_values, NTICKS)

        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            SIR.Susceptible(model),
            SIR.Infectious(model, infdurdist),
            SIR.Recovered(model),
            TransmissionSI(model, infdurdist, seasonality=seasonality),
        ]

        model.run("SIR Spatial Seasonality")

        # Calculate peak times for each node
        peak_times = []
        for node_idx in range(3):
            I_series = model.nodes.I[:, node_idx]
            peak_time = np.argmax(I_series)
            peak_times.append(peak_time)

        attenuated_peak, baseline_peak, amplified_peak = peak_times

        # With effective R0 = 0.875 (0.5 * 0.25 * 7), attenuated node should have early, low peak
        # Baseline R0 = 1.75, amplified R0 = 3.5 should have earlier peaks due to faster spread
        assert amplified_peak < baseline_peak, (
            f"Amplified node peak (day {amplified_peak}) should occur before baseline peak (day {baseline_peak})"
        )

        # Calculate final attack rates
        attack_rates = []
        for node_idx in range(3):
            attack_rate = model.nodes.R[-1, node_idx] / POPULATION
            attack_rates.append(attack_rate)

        attenuated_ar, baseline_ar, amplified_ar = attack_rates

        # Amplified should have higher attack rate than baseline
        assert amplified_ar > baseline_ar, f"Amplified node attack rate ({amplified_ar:.4f}) should exceed baseline ({baseline_ar:.4f})"

        # Attenuated with R0 < 1 should have much lower attack rate
        assert attenuated_ar < baseline_ar, (
            f"Attenuated node attack rate ({attenuated_ar:.4f}) should be less than baseline ({baseline_ar:.4f})"
        )

        return

    def test_seir_spatial_seasonality(self):
        """
        Validate spatially independent seasonal forcing across disconnected nodes (SEIR model).

        WHAT IS VALIDATED:
        - E peaks before I for nodes with R₀>1 (baseline and amplified nodes)
        - Amplified node peaks before baseline node
        - Attack rates follow ordering: amplified > baseline > attenuated
        - Attenuated node (R₀<1) has minimal attack rate

        SCENARIO SETUP:
        - Three nodes (1×3 grid), each with 100,000 population, 500 initial infections
        - Beta = 0.25, exp=5d, inf=7d, seasonality = [0.5×, 1.0×, 2.0×]
        - Effective R₀ = [0.875, 1.75, 3.5] across nodes
        - gravity_k = 0.0 (nodes isolated)
        - Skips E-before-I check for node 0 (immediate die-out)

        FAILURE MEANING:
        If this test fails, spatial seasonality isn't correctly applied to SEIR transmission. Each
        node should exhibit independent SEIR dynamics with its assigned R₀. Node 0 epidemic dies out
        immediately (R₀<1), so its I peak is at t=0 - checking E-before-I there would spuriously fail.
        Failure in E-before-I for nodes 1-2 indicates broken SEIR progression. Wrong attack rate or
        timing order suggests seasonality not applied per-node or values mapped incorrectly.
        """

        # Create 3-node scenario with no spatial transmission
        scenario = stdgrid(M=1, N=3, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 500
        scenario["E"] = 0
        scenario["I"] = 500
        scenario["R"] = 0

        beta = 0.25
        exp_mean = 5.0
        inf_mean = 7.0

        # Disable spatial transmission by setting gravity_k to 0
        params = PropertySet({"nticks": NTICKS, "beta": beta, "gravity_k": 0.0})
        model = Model(scenario, params)

        # Create spatially varying seasonality: [0.5, 1.0, 2.0]
        seasonality_values = np.array([0.5, 1.0, 2.0], dtype=np.float32)
        seasonality = ValuesMap.from_nodes(seasonality_values, NTICKS)

        expdurdist = dists.normal(loc=exp_mean, scale=1.0)
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            SEIR.Susceptible(model),
            SEIR.Exposed(model, expdurdist, infdurdist),
            SEIR.Infectious(model, infdurdist),
            SEIR.Recovered(model),
            TransmissionSE(model, expdurdist, seasonality=seasonality),
        ]

        model.run("SEIR Spatial Seasonality")

        # For nodes with R0 > 1 (baseline and amplified), verify E peaks before I
        # Skip attenuated node (node 0) where epidemic dies out immediately
        for node_idx in [1, 2]:  # Baseline and amplified nodes only
            E_series = model.nodes.E[:, node_idx]
            I_series = model.nodes.I[:, node_idx]
            E_peak_time = np.argmax(E_series)
            I_peak_time = np.argmax(I_series)

            assert E_peak_time < I_peak_time, f"Node {node_idx}: E should peak (day {E_peak_time}) before I (day {I_peak_time})"

        # Calculate I peak times for each node
        I_peak_times = []
        for node_idx in range(3):
            I_series = model.nodes.I[:, node_idx]
            peak_time = np.argmax(I_series)
            I_peak_times.append(peak_time)

        attenuated_peak, baseline_peak, amplified_peak = I_peak_times

        # Amplified node should peak earlier than baseline
        assert amplified_peak < baseline_peak, (
            f"Amplified node I peak (day {amplified_peak}) should occur before baseline peak (day {baseline_peak})"
        )

        # Calculate final attack rates
        attack_rates = []
        for node_idx in range(3):
            attack_rate = model.nodes.R[-1, node_idx] / POPULATION
            attack_rates.append(attack_rate)

        attenuated_ar, baseline_ar, amplified_ar = attack_rates

        # Amplified should have higher attack rate than baseline
        assert amplified_ar > baseline_ar, f"Amplified node attack rate ({amplified_ar:.4f}) should exceed baseline ({baseline_ar:.4f})"

        # Attenuated with R0 < 1 should have much lower attack rate
        assert attenuated_ar < baseline_ar, (
            f"Attenuated node attack rate ({attenuated_ar:.4f}) should be less than baseline ({baseline_ar:.4f})"
        )

        return

    # ========================================================================
    # Temporally Varying Seasonality Tests
    # ========================================================================

    def test_si_temporal_seasonality(self):
        """
        Validate time-varying seasonal forcing that declines during outbreak (SI model).

        WHAT IS VALIDATED:
        - Final infections are lower with declining seasonality vs baseline
        - Growth during decline period (25%-75% of infections) is reduced vs baseline
        - Declining transmission from 1.0→0.0 slows epidemic progression

        SCENARIO SETUP:
        - Two-phase approach:
            (1) baseline run to find 25%/75% infection thresholds
            (2) seasonal run with seasonality declining from 1.0 to 0.0 between those times
        - Beta = 0.2 (R₀ between 1-2 for slower outbreak), 1000-day simulation
        - Seasonality = 1.0 until 25%, linearly declines to 0.0 by 75%, stays 0.0 after
        - 100 initial infections (small seed for gradual spread)

        FAILURE MEANING:
        If this test fails, temporal seasonality isn't modulating transmission over time. The
        declining seasonality should slow epidemic growth during the 25%-75% window, reducing
        final infections compared to constant transmission. Failure indicates seasonality not
        being updated each timestep, ValuesMap.from_timeseries() not working, or transmission
        calculation not reading current seasonality value. Equal growth suggests seasonality ignored.
        """

        # Use lower beta and longer simulation to see temporal effects
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 100
        scenario["I"] = 100

        # Lower R0 (between 1 and 2) for slower outbreak
        beta = 0.2  # With SI model, this gives gradual spread
        long_nticks = 1000  # Longer simulation to see full outbreak

        params = PropertySet({"nticks": long_nticks, "beta": beta})
        model = Model(scenario, params)

        # Start with seasonality = 1.0, will update after first run
        temp_seasonality = ValuesMap.from_scalar(1.0, long_nticks, 1)

        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=temp_seasonality),
        ]

        # First run to determine when cumulative infections reach 25% and 75%
        model.run("SI Temporal - Initial Run")

        # Calculate cumulative new infections (change in I + change in recovered if we had it)
        # For SI model, cumulative infections = current I count
        cumulative_infections = model.nodes.I[:, 0]
        total_infections = cumulative_infections[-1]

        # Find when we reach 25% and 75% of total infections
        threshold_25 = total_infections * 0.25
        threshold_75 = total_infections * 0.75

        time_25 = np.argmax(cumulative_infections >= threshold_25)
        time_75 = np.argmax(cumulative_infections >= threshold_75)

        # Create declining seasonality between these times
        seasonality_array = np.ones(long_nticks, dtype=np.float32)
        decline_length = time_75 - time_25
        if decline_length > 0:
            seasonality_array[time_25:time_75] = np.linspace(1.0, 0.0, decline_length)
            seasonality_array[time_75:] = 0.0

        seasonality = ValuesMap.from_timeseries(seasonality_array, 1)

        # Second run with declining seasonality
        model2 = Model(scenario, params)
        model2.components = [
            SI.Susceptible(model2),
            SI.Infectious(model2),
            TransmissionSIx(model2, seasonality=seasonality),
        ]

        model2.run("SI Temporal - With Declining Seasonality")

        # Verify that declining seasonality slowed the outbreak
        final_I_baseline = model.nodes.I[-1, 0]
        final_I_seasonal = model2.nodes.I[-1, 0]

        # With declining transmission, final infections should be lower
        assert final_I_seasonal < final_I_baseline, (
            f"Declining seasonality should reduce final infections: seasonal={final_I_seasonal}, baseline={final_I_baseline}"
        )

        # Verify that the decline period shows reduced growth
        growth_baseline = model.nodes.I[time_75, 0] - model.nodes.I[time_25, 0]
        growth_seasonal = model2.nodes.I[time_75, 0] - model2.nodes.I[time_25, 0]

        assert growth_seasonal < growth_baseline, (
            f"Growth during decline period should be reduced: seasonal={growth_seasonal}, baseline={growth_baseline}"
        )

        return

    def test_sir_temporal_seasonality(self):
        """
        Validate time-varying seasonal forcing that declines during outbreak (SIR model).

        WHAT IS VALIDATED:
        - Final attack rate is lower with declining seasonality vs baseline
        - Cumulative infections at 75% threshold are lower than baseline
        - Declining seasonality interrupts natural SIR epidemic progression

        SCENARIO SETUP:
        - Two-phase:
            (1) baseline run to find 25%/75% cumulative infection thresholds
            (2) seasonal run with seasonality 1.0→0.0 between those times
        - Beta = 0.2, inf_period = 7d (R₀ ≈ 1.4), 1500-day simulation
        - Cumulative infections = I + R (total ever infected)
        - 100 initial infections for gradual spread

        FAILURE MEANING:
        If this test fails, temporal seasonality isn't properly reducing SIR transmission over time.
        The declining force should reduce the final attack rate below what natural SIR dynamics
        would achieve. Failure indicates time-varying seasonality not applied each timestep, or the
        decline isn't strong enough to overcome existing momentum. Equal attack rates suggest
        seasonality parameter not being read or applied in TransmissionSI step function.
        """

        # Use R0 between 1 and 2 for slower outbreak
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 100
        scenario["I"] = 100
        scenario["R"] = 0

        beta = 0.2  # Lower beta for R0 ~ 1.4 (0.2 * 7 = 1.4)
        inf_mean = 7.0
        long_nticks = 1500

        params = PropertySet({"nticks": long_nticks, "beta": beta})
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        # First run without seasonality to determine timing
        model = Model(scenario, params)
        model.components = [
            SIR.Susceptible(model),
            SIR.Infectious(model, infdurdist),
            SIR.Recovered(model),
            TransmissionSI(model, infdurdist, seasonality=None),
        ]

        model.run("SIR Temporal - Initial Run")

        # Calculate cumulative infections (sum of I + R over time)
        cumulative_infections = model.nodes.I[:, 0] + model.nodes.R[:, 0]
        total_infections = cumulative_infections[-1]

        threshold_25 = total_infections * 0.25
        threshold_75 = total_infections * 0.75

        time_25 = np.argmax(cumulative_infections >= threshold_25)
        time_75 = np.argmax(cumulative_infections >= threshold_75)

        # Create declining seasonality
        seasonality_array = np.ones(long_nticks, dtype=np.float32)
        decline_length = time_75 - time_25
        if decline_length > 0:
            seasonality_array[time_25:time_75] = np.linspace(1.0, 0.0, decline_length)
            seasonality_array[time_75:] = 0.0

        seasonality = ValuesMap.from_timeseries(seasonality_array, 1)

        # Second run with declining seasonality
        model2 = Model(scenario, params)
        model2.components = [
            SIR.Susceptible(model2),
            SIR.Infectious(model2, infdurdist),
            SIR.Recovered(model2),
            TransmissionSI(model2, infdurdist, seasonality=seasonality),
        ]

        model2.run("SIR Temporal - With Declining Seasonality")

        # Verify that declining seasonality reduced final attack rate
        final_ar_baseline = model.nodes.R[-1, 0] / POPULATION
        final_ar_seasonal = model2.nodes.R[-1, 0] / POPULATION

        assert final_ar_seasonal < final_ar_baseline, (
            f"Declining seasonality should reduce final attack rate: seasonal={final_ar_seasonal:.4f}, baseline={final_ar_baseline:.4f}"
        )

        # Check cumulative infections at time_75
        cum_inf_baseline_75 = cumulative_infections[time_75]
        cum_inf_seasonal_75 = model2.nodes.I[time_75, 0] + model2.nodes.R[time_75, 0]

        assert cum_inf_seasonal_75 < cum_inf_baseline_75, (
            f"Cumulative infections at 75% mark should be lower with declining seasonality: "
            f"seasonal={cum_inf_seasonal_75}, baseline={cum_inf_baseline_75}"
        )

        return

    def test_seir_temporal_seasonality(self):
        """
        Validate time-varying seasonal forcing that declines during outbreak (SEIR model).

        WHAT IS VALIDATED:
        - Final attack rate is lower with declining seasonality vs baseline
        - E peaks before I (if epidemic develops beyond initial seed)
        - Cumulative infections at 75% threshold are lower than baseline
        - Declining seasonality interrupts natural SEIR progression

        SCENARIO SETUP:
        - Two-phase:
            (1) baseline run to find 25%/75% cumulative infection thresholds
            (2) seasonal run with seasonality 1.0→0.0 between those times
        - Beta = 0.2, exp=5d, inf=7d (R₀ ≈ 1.4), 1500-day simulation
        - Cumulative infections = E + I + R (total who left S)
        - 100 initial infections for gradual spread
        - E-before-I check only if I_peak > 10 (epidemic actually develops)

        FAILURE MEANING:
        If this test fails, temporal seasonality isn't properly reducing SEIR transmission over time.
        The declining force should slow E→I→R progression and reduce final attack rate. Failure
        indicates time-varying seasonality not applied in TransmissionSE, or compartment dynamics
        don't respond to changing transmission. Equal attack rates mean seasonality ignored. E-before-I
        failure (when I_peak>10) indicates broken SEIR compartment progression independent of seasonality.
        """

        # Use R0 between 1 and 2 for slower outbreak
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 100
        scenario["E"] = 0
        scenario["I"] = 100
        scenario["R"] = 0

        beta = 0.2  # Lower beta for R0 ~ 1.4 (0.2 * 7 = 1.4)
        exp_mean = 5.0
        inf_mean = 7.0
        long_nticks = 1500

        params = PropertySet({"nticks": long_nticks, "beta": beta})
        expdurdist = dists.normal(loc=exp_mean, scale=1.0)
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        # First run without seasonality to determine timing
        model = Model(scenario, params)
        model.components = [
            SEIR.Susceptible(model),
            SEIR.Exposed(model, expdurdist, infdurdist),
            SEIR.Infectious(model, infdurdist),
            SEIR.Recovered(model),
            TransmissionSE(model, expdurdist, seasonality=None),
        ]

        model.run("SEIR Temporal - Initial Run")

        # Calculate cumulative infections (E + I + R over time)
        cumulative_infections = model.nodes.E[:, 0] + model.nodes.I[:, 0] + model.nodes.R[:, 0]
        total_infections = cumulative_infections[-1]

        threshold_25 = total_infections * 0.25
        threshold_75 = total_infections * 0.75

        time_25 = np.argmax(cumulative_infections >= threshold_25)
        time_75 = np.argmax(cumulative_infections >= threshold_75)

        # Create declining seasonality
        seasonality_array = np.ones(long_nticks, dtype=np.float32)
        decline_length = time_75 - time_25
        if decline_length > 0:
            seasonality_array[time_25:time_75] = np.linspace(1.0, 0.0, decline_length)
            seasonality_array[time_75:] = 0.0

        seasonality = ValuesMap.from_timeseries(seasonality_array, 1)

        # Second run with declining seasonality
        model2 = Model(scenario, params)
        model2.components = [
            SEIR.Susceptible(model2),
            SEIR.Exposed(model2, expdurdist, infdurdist),
            SEIR.Infectious(model2, infdurdist),
            SEIR.Recovered(model2),
            TransmissionSE(model2, expdurdist, seasonality=seasonality),
        ]

        model2.run("SEIR Temporal - With Declining Seasonality")

        # Verify that declining seasonality reduced final attack rate
        final_ar_baseline = model.nodes.R[-1, 0] / POPULATION
        final_ar_seasonal = model2.nodes.R[-1, 0] / POPULATION

        assert final_ar_seasonal < final_ar_baseline, (
            f"Declining seasonality should reduce final attack rate: seasonal={final_ar_seasonal:.4f}, baseline={final_ar_baseline:.4f}"
        )

        # Verify E peaks before I in the seasonal model (if epidemic actually grows)
        E_peak_seasonal = np.argmax(model2.nodes.E[:, 0])
        I_peak_seasonal = np.argmax(model2.nodes.I[:, 0])

        # Only check if I_peak is not at the start (epidemic actually develops)
        if I_peak_seasonal > 10:
            assert E_peak_seasonal < I_peak_seasonal, (
                f"E should peak before I even with declining seasonality: E_peak={E_peak_seasonal}, I_peak={I_peak_seasonal}"
            )

        # Verify reduced cumulative infections at time_75
        cum_inf_baseline_75 = cumulative_infections[time_75]
        cum_inf_seasonal_75 = model2.nodes.E[time_75, 0] + model2.nodes.I[time_75, 0] + model2.nodes.R[time_75, 0]

        assert cum_inf_seasonal_75 < cum_inf_baseline_75, (
            f"Cumulative infections at 75% mark should be lower with declining seasonality: "
            f"seasonal={cum_inf_seasonal_75}, baseline={cum_inf_baseline_75}"
        )

        return

    # ========================================================================
    # Helper Methods
    # ========================================================================

    def _create_si_baseline(self):
        """Create baseline SI model for comparison."""
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 1000
        scenario["I"] = 1000

        params = PropertySet({"nticks": NTICKS, "beta": 0.5})
        model = Model(scenario, params)

        model.components = [
            SI.Susceptible(model),
            SI.Infectious(model),
            TransmissionSIx(model, seasonality=None),
        ]

        model.run("SI Baseline")
        return model


def plot_models(*models):
    """Utility to plot model compartments for debugging."""
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 9))

    traces = [("S", "Susceptible", "blue"), ("E", "Exposed", "orange"), ("I", "Infectious", "red"), ("R", "Recovered", "green")]
    styles = ["-", "--", "-.", ":"]

    for i, model in enumerate(models):
        for comp, label, color in traces:
            if hasattr(model.nodes, comp):
                plt.plot(
                    getattr(model.nodes, comp).sum(axis=1),
                    label=f"{label} (Model {i + 1})",
                    color=color,
                    linestyle=styles[i % len(styles)],
                )

    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.title("Model Compartments Over Time")
    plt.legend()
    plt.grid()
    plt.show()

    return


if __name__ == "__main__":
    unittest.main()
