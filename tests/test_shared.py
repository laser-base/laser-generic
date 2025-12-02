import unittest

import numpy as np
from laser.core.demographics import AliasedDistribution
from laser.core.demographics import KaplanMeierEstimator
from scipy.stats import ks_2samp

from laser.generic.shared import sample_dobs
from laser.generic.shared import sample_dods


class TestShared(unittest.TestCase):
    def test_sample_dobs_with_tick_zero_uniform_pyramid(self):
        pyramid = AliasedDistribution([1000] * 100)  # Uniform distribution over ages 0..99

        n_agents = 100_000
        dobs = np.zeros(n_agents, dtype=np.int32)
        tick = 0

        sample_dobs(pyramid, dobs, tick)

        # The minimum dob should be at least -99*365-364 (oldest, max noise)
        assert np.all(dobs > -100 * 365), "Minimum dob is less than expected."
        # The maximum dob should be at most 0 (youngest, min noise)
        assert np.all(dobs <= 0), "Maximum dob is greater than expected."

        # Perform a KS test to verify that the dob distribution matches the population pyramid

        # Convert dobs to ages in years (since dobs are negative ages in days)
        ages_sampled = -dobs // 365

        # Generate the expected population ages according to the pyramid
        expected_ages = np.repeat(np.arange(100), np.round(n_agents * pyramid.probs / pyramid.probs.sum()).astype(np.int32))

        # Perform a chi-squared test to verify that the dob distribution matches the population pyramid
        observed_counts = np.bincount(ages_sampled, minlength=100)
        expected_counts = np.bincount(expected_ages, minlength=100)
        from scipy.stats import chisquare

        chi2_stat, p_value = chisquare(observed_counts, expected_counts)

        assert p_value > 0.01, f"Chi-squared test failed: p-value={p_value}"

        return

    def test_sample_dobs_with_tick_zero_exponential_pyramid(self):
        # Exponential decay: P(age) ~ exp(-lambda * age)
        n_bins = 100
        n_agents = 100_000
        decay_rate = 0.05  # Chosen so that population decreases noticeably over 100 bins

        # Create exponential population pyramid
        ages = np.arange(n_bins)
        weights = np.exp(-decay_rate * ages)
        weights /= weights.sum()  # Normalize to probabilities
        counts = np.round(weights * n_agents)
        pyramid = AliasedDistribution(counts.astype(np.int32))

        dobs = np.zeros(n_agents, dtype=np.int32)
        tick = 0

        sample_dobs(pyramid, dobs, tick)

        # The minimum dob should be at least -99*365-364 (oldest, max noise)
        assert np.all(dobs > -100 * 365), "Minimum dob is less than expected."
        # The maximum dob should be at most 0 (youngest, min noise)
        assert np.all(dobs <= 0), "Maximum dob is greater than expected."

        # Convert dobs to ages in years
        ages_sampled = -dobs // 365

        # Generate the expected population ages according to the pyramid
        expected_ages = np.repeat(np.arange(n_bins), (weights * n_agents).astype(int))

        # Perform the two-sample KS test
        _ks_stat, p_value = ks_2samp(ages_sampled, expected_ages)

        assert p_value > 0.01, f"KS test failed: p-value={p_value}"

        return

    def test_sample_dods_with_uniform_hazard(self):
        n_agents = 100_000
        dobs = np.zeros(n_agents, dtype=np.int32)
        tick = 0

        # Create a KaplanMeierEstimator with equal hazard: 1000 deaths per year for 100 years
        n_years = 100
        deaths_per_year = 1000
        hazards = np.full(n_years, deaths_per_year)
        estimator = KaplanMeierEstimator(hazards.cumsum())

        dods = np.zeros(n_agents, dtype=np.int32)
        sample_dods(dobs, estimator, tick, dods)

        # The maximum dod should be at most 99*365+364 (oldest, max noise)
        assert np.all(dods <= 100 * 365), "Maximum dod is greater than expected."
        # The minimum dod should be at least 0 (youngest, min noise)
        assert np.all(dods >= 0), "Minimum dod is less than expected."

        # Convert dods to ages at death in years
        ages_at_death = dods // 365

        # Generate the expected death ages according to the estimator
        # For uniform hazard, the probability of dying at each year is equal
        expected_ages = np.repeat(np.arange(n_years), n_agents // n_years)

        # Perform the two-sample KS test
        _ks_stat, p_value = ks_2samp(ages_at_death, expected_ages)

        assert p_value > 0.01, f"KS test failed: p-value={p_value}"

        return

    def test_sample_dods_survival_curve_matches_estimator(self):
        n_agents = 100_000
        n_years = 100
        deaths_per_year = 1000
        hazards = np.full(n_years, deaths_per_year)
        estimator = KaplanMeierEstimator(hazards.cumsum())

        # Uniform population pyramid over ages 0..99
        pyramid = AliasedDistribution([1_000] * n_years)
        dobs = np.zeros(n_agents, dtype=np.int32)
        tick = 0
        sample_dobs(pyramid, dobs, tick)

        dods = np.zeros(n_agents, dtype=np.int32)
        sample_dods(dobs, estimator, tick, dods)

        # For each age group, compare empirical survival curve to estimator's survival curve
        ages = -dobs // 365
        ages_at_death = (dods - dobs) // 365

        for age in range(0, n_years, 10):  # Check every 10 years for speed
            mask = ages == age
            if np.sum(mask) < 100:
                continue  # Not enough samples for this age group

            counts = np.bincount(ages_at_death[mask])
            # Counts should range from age to n_years-1
            empirical = counts.cumsum() / counts.sum()

            offset = estimator.cumulative_deaths[age - 1] if age > 0 else 0
            expected = (estimator.cumulative_deaths[age:n_years] - offset) / (estimator.cumulative_deaths[-1] - offset)

            empirical = empirical[-len(expected) :]  # Align lengths

            # Compare curves: mean absolute error should be small
            mae = np.mean(np.abs(np.array(empirical) - expected))
            assert mae < 0.05, f"Survival curve mismatch at age {age}: MAE={mae}"

        return


if __name__ == "__main__":
    unittest.main()
