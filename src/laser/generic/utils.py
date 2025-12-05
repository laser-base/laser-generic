"""
This module provides utility functions for the laser-measles project.

"""

from math import ceil

from laser.core import PropertySet


def get_default_parameters() -> PropertySet:
    """
    Returns a default PropertySet with common parameters used across laser-generic models.

    Each parameter in the returned PropertySet is described below, along with its default value and rationale:

        nticks (int, default=730): Number of simulation ticks (days). Default is 2 years (365*2), which is a typical duration for seasonal epidemic simulations.
        beta (float, default=0.15): Transmission rate per contact. Chosen as a moderate value for SIR-type models to reflect realistic disease spread.
        biweekly_beta_scalar (list of float, default=[1.0]*biweekly_steps): Scalar for beta for each biweekly period. Default is 1.0 for all periods, meaning no seasonal variation unless specified.
        cbr (float, default=0.03): Constant birth rate. Set to 0.03 to represent a typical annual birth rate in population models.
        exp_shape (float, default=2.0): Shape parameter for the exposed period distribution. Default chosen for moderate dispersion.
        exp_scale (float, default=2.0): Scale parameter for the exposed period distribution. Default chosen for moderate mean duration.
        inf_mean (float, default=4.0): Mean infectious period (days). Set to 4.0 to reflect typical infectious durations for diseases like measles.
        inf_sigma (float, default=1.0): Standard deviation of infectious period. Default is 1.0 for moderate variability.
        seasonality_factor (float, default=0.2): Amplitude of seasonal forcing. Chosen to allow moderate seasonal variation in transmission.
        seasonality_phase (float, default=0.0): Phase offset for seasonality. Default is 0.0, meaning no phase shift.
        importation_count (int, default=1): Number of cases imported per importation event. Default is 1 for sporadic importation.
        importation_period (int, default=30): Days between importation events. Default is 30 to represent monthly importation.
        importation_start (int, default=0): Start day for importation events. Default is 0 (simulation start).
        importation_end (int, default=730): End day for importation events. Default is 2 years (365*2).
        seed (int, default=123): Random seed for reproducibility. Default is 123.
        verbose (bool, default=False): If True, enables verbose output. Default is False for minimal output.
    These values are chosen to be broadly reasonable for seasonal SIR-type models with importation.

    We need a function like this because even-though laser-core requires no particular param name,
    laser-generic code does presume certain parameters and there's no elegant way to just discover
    what those are. So we put them here.
    """
    nticks = 365 * 2
    biweekly_steps = ceil(nticks / 14)
    return PropertySet(
        {
            "nticks": nticks,
            "beta": 0.15,
            "biweekly_beta_scalar": [1.0] * biweekly_steps,
            "cbr": 0.03,
            "exp_shape": 2.0,
            "exp_scale": 2.0,
            "inf_mean": 4.0,
            "inf_sigma": 1.0,
            "seasonality_factor": 0.2,
            "seasonality_phase": 0.0,
            "importation_count": 1,
            "importation_period": 30,
            "importation_start": 0,
            "importation_end": 365 * 2,
            "seed": 123,
            "verbose": False,
        }
    )
