import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from laser.core import PropertySet
from laser.generic.model import Model
from laser.generic.utils import get_default_parameters
from laser.generic.newutils import grid
from laser.generic.components import (
    Susceptible,
    InfectiousSI,
    TransmissionSI,
    State,
)

def si_logistic(t, beta, N, t0):
    """Analytical SI logistic curve for total infections."""
    return N / (1 + (N - 1) * np.exp(-beta * (t - t0)))


def test_si_logistic_beta_recovery():
    """
    Logistic-fitting test for pure SI model using the documented SI module.

    This matches the real LASER SI design:
    - No recovery
    - No itimer required
    - S -> I transitions only
    """

    from laser.generic import SI

    pop = 100_000
    n_days = 200
    betas = [0.01, 0.02, 0.03]

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    output = []

    for beta in betas:
        params = get_default_parameters() | {
            "nticks": n_days,
            "seed": 42,
            "beta": beta,
        }

        # Initialize scenario S/I explicitly
        scenario["S"] = pop - 1
        scenario["I"] = 1

        model = Model(scenario, params)

        # Correct SI components — NO itimer needed
        S = SI.Susceptible(model)
        I = SI.Infectious(model)
        T = SI.Transmission(model)

        model.components = [S, I, T]

        model.run("SI logistic test")

        I_series = model.nodes.I[:, 0]
        t = np.arange(len(I_series))

        popt, _ = curve_fit(
            si_logistic,
            t,
            I_series,
            p0=[beta * 0.9, pop, 5],
            maxfev=5000,
        )

        fitted_beta = popt[0]
        output.append((beta, fitted_beta))

    # Assert β recovery within 5%
    for true_beta, fitted_beta in output:
        rel_error = abs(fitted_beta - true_beta) / true_beta
        assert rel_error < 0.05, \
            f"β recovery failed: true={true_beta}, fitted={fitted_beta}"
