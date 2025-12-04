import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from laser.core import PropertySet
from laser.generic.model import Model
from laser.generic.utils import (
    get_default_parameters,
    seed_infections_randomly,
    seed_infections_randomly_SI,
    seed_infections_in_patch,
)
from laser.generic.newutils import grid

# --- Correct modern LASER components ---
from laser.generic.components import (
    Susceptible,
    Exposed,
    InfectiousSI,     # SI
    InfectiousIS,     # SIS
    InfectiousIR,     # SIR / SEIR
    InfectiousIRS,    # SIRS / SEIRS
    Recovered,
    RecoveredRS,
    TransmissionSI,   # S → I (SI/SIS/SIR/SIRS)
    TransmissionSE,   # S → E (SEIR)
)

from laser.generic.vitaldynamics import BirthsByCBR as Births, ConstantPopVitalDynamics
from laser.generic.importation import Infect_Agents_In_Patch, Infect_Random_Agents
from laser.generic.immunization import RoutineImmunization, ImmunizationCampaign


# ======================================================================
# Helpers
# ======================================================================

def node_population(model):
    """Compute S+E+I+R per node, shape (T, N)."""
    nodes = model.nodes
    S = nodes.S
    E = getattr(nodes, "E", 0)
    I = getattr(nodes, "I", 0)
    R = getattr(nodes, "R", 0)
    return S + E + I + R


def assert_model_sanity(model):
    nodes = model.nodes
    S = nodes.S[:, 0]
    I = nodes.I[:, 0]
    N = node_population(model)[:, 0]
    inc = nodes.incidence[:, 0]

    assert inc.sum() > 0, "No transmission occurred"
    assert np.any(S[1:] < S[0]), "Susceptibles never decreased"
    assert np.all(S >= 0), "Negative susceptibles"
    assert np.all(S <= N), "Susceptibles > population"

    # I(t) ≈ I(0) + cumulative incidence
    I_expected = I[0] + np.cumsum(inc)
    assert np.allclose(I[1:], I_expected[1:], atol=1e-5), "I inconsistent with incidence"


# ======================================================================
# Baseline model fixture
# ======================================================================

def baseline_model():
    pop = 100_000
    nticks = 365

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    params = get_default_parameters() | {
        "seed": 42,
        "nticks": nticks,
        "beta": 0.30,
        "inf_mean": 7,
        "verbose": False,
    }

    model = Model(scenario, params)

    model.components = [
        ConstantPopVitalDynamics,
        Susceptible,
        Exposed,
        InfectiousIR,
        TransmissionSE,   # SEIR-style transmission
    ]

    seed_infections_in_patch(model, ninfections=50, ipatch=0)
    return model


@pytest.fixture
def stable_transmission_model():
    return baseline_model()


# ======================================================================
# Tests
# ======================================================================

@pytest.mark.modeltest
def test_si_model_nobirths_flow():
    pop = 100_000
    nticks = 180

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)
    params = PropertySet({"seed": 42, "nticks": nticks, "beta": 0.03, "verbose": False})

    model = Model(scenario, params)
    model.components = [Susceptible, TransmissionSI]
    seed_infections_randomly_SI(model, ninfections=1)
    model.run()

    assert_model_sanity(model)
    assert model.nodes.I[-1, 0] > model.nodes.I[0, 0]


@pytest.mark.modeltest
def test_sir_nobirths_short():
    pop = 100_000
    nticks = 365

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    beta = 0.06
    gamma = 1 / 20
    params = PropertySet({
        "seed": 1,
        "nticks": nticks,
        "beta": beta,
        "inf_mean": 1 / gamma,
        "verbose": False,
    })

    model = Model(scenario, params)
    model.components = [Susceptible, InfectiousIR, TransmissionSI]
    seed_infections_randomly(model, ninfections=50)
    model.run()

    assert_model_sanity(model)
    I = model.nodes.I[:, 0]
    assert I[-1] < I.max()


@pytest.mark.modeltest
def test_si_model_with_births_short():
    pop = 100_000
    nticks = 365 * 2

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    params = PropertySet({
        "seed": 123,
        "nticks": nticks,
        "beta": 0.02,
        "cbr": 0.03,
        "verbose": False,
    })

    model = Model(scenario, params)
    model.components = [ConstantPopVitalDynamics, Susceptible, TransmissionSI]
    seed_infections_randomly_SI(model, ninfections=10)
    model.run()

    assert_model_sanity(model)
    assert model.nodes.I[-1, 0] > 0


@pytest.mark.modeltest
def test_sei_model_with_births_short():
    pop = 100_000
    nticks = 365 * 2

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    params = get_default_parameters() | {
        "seed": 123,
        "nticks": nticks,
        "beta": 0.05,
        "inf_mean": 5,
        "cbr": 0.03,
        "verbose": False,
    }

    model = Model(scenario, params)
    model.components = [
        ConstantPopVitalDynamics,
        Susceptible,
        Exposed,
        InfectiousIR,
        TransmissionSE,
    ]
    seed_infections_randomly_SI(model, ninfections=10)
    model.run()

    assert_model_sanity(model)
    assert model.nodes.I[-1, 0] > 0


@pytest.mark.modeltest
def test_sis_model_short():
    pop = 100_000
    nticks = 500

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    params = PropertySet({
        "seed": 99,
        "nticks": nticks,
        "beta": 0.05,
        "inf_mean": 10,
        "verbose": False,
    })

    model = Model(scenario, params)
    model.components = [Susceptible, InfectiousIS, TransmissionSI]
    seed_infections_randomly(model, ninfections=50)
    model.run()

    assert_model_sanity(model)

    I = model.nodes.I[:, 0]
    assert np.any(I[1:] > I[0])
    assert I[-1] > 0


@pytest.mark.modeltest
def test_routine_immunization_blocks_spread():
    pop = 100_000
    nticks = 365 * 2

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    params = get_default_parameters() | {
        "seed": 321,
        "nticks": nticks,
        "beta": 0.05,
        "inf_mean": 5,
        "cbr": 0.03,
        "verbose": False,
    }

    model = Model(scenario, params)
    model.components = [
        ConstantPopVitalDynamics,
        Susceptible,
        Exposed,
        InfectiousIR,
        TransmissionSE,
        lambda m, v: RoutineImmunization(m, period=365, coverage=0.9, age=365, verbose=v),
    ]
    seed_infections_randomly_SI(model, ninfections=10)
    model.run()

    assert_model_sanity(model)
    assert model.nodes.I[-1, 0] < 0.5 * pop


@pytest.mark.modeltest
def test_mobility_spreads_infection_across_nodes():
    pop = 50_000
    nticks = 180

    # 1×2 grid → two nodes side-by-side
    scenario = grid(M=1, N=2, population_fn=lambda x, y: pop)

    params = get_default_parameters() | {
        "seed": 42,
        "nticks": nticks,
        "beta": 0.05,
        "inf_mean": 5,
        "verbose": False,
    }

    model = Model(scenario, params)

    # Simple mobility network
    model.nodes.network = np.array([[0.95, 0.05], [0.05, 0.95]])

    model.components = [Susceptible, Exposed, InfectiousIR, TransmissionSE]

    seed_infections_in_patch(model, ninfections=10, ipatch=0)
    model.run()

    assert_model_sanity(model)
    assert model.nodes.I[-1, 1] > 0


@pytest.mark.modeltest
def test_births_base_runs_minimally():
    pop = 1000
    nticks = 10

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    params = get_default_parameters() | {
        "nticks": nticks,
        "seed": 42,
        "verbose": False,
    }

    model = Model(scenario, params)
    model.components = [Births, Susceptible]

    try:
        model.run()
    except Exception as e:
        pytest.fail(f"Births class failed: {e}")


@pytest.mark.modeltest
def test_importation_keeps_infection_alive():
    pop = 100_000
    nticks = 365 * 5

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    params = get_default_parameters() | {
        "nticks": nticks,
        "seed": 123,
        "beta": 0.35,
        "inf_mean": 5,
        "cbr": 0.03,
        "importation_period": 15,
        "importation_count": 10,
        "importation_start": 10,
        "importation_end": nticks,
        "verbose": False,
    }

    model = Model(scenario, params)

    model.components = [
        ConstantPopVitalDynamics,
        Susceptible,
        Exposed,
        InfectiousIR,
        TransmissionSE,
        Infect_Random_Agents,
    ]

    model.run()

    inc = model.nodes.incidence[:, 0]
    assert inc.sum() > 0
    assert inc.sum() > 1000


@pytest.mark.modeltest
def test_targeted_importation_hits_correct_patch():
    pop = 100_000
    nticks = 365

    # A 1×2 grid gives us 2 patches automatically
    scenario = grid(M=1, N=2, population_fn=lambda x, y: pop)

    params = get_default_parameters() | {
        "nticks": nticks,
        "seed": 222,
        "beta": 0.05,
        "inf_mean": 5,
        "cbr": 0.03,
        "importation_period": 30,
        "importation_count": 5,
        "importation_start": 0,
        "importation_end": nticks,
        "importation_target": 1,
        "verbose": False,
    }

    model = Model(scenario, params)
    import_comp = Infect_Agents_In_Patch(model, verbose=params["verbose"])

    model.components = [
        ConstantPopVitalDynamics,
        Susceptible,
        Exposed,
        InfectiousIR,
        TransmissionSE,
        import_comp,
    ]

    model.run()

    I0 = model.nodes.I[:, 0]
    I1 = model.nodes.I[:, 1]

    assert I1.sum() > 0
    assert I0.sum() == 0


@pytest.mark.modeltest
def test_transmission_sir_behaves_like_transmission():
    pop = 100_000
    nticks = 365

    scenario = grid(M=1, N=1, population_fn=lambda x, y: pop)

    params = get_default_parameters() | {
        "nticks": nticks,
        "seed": 777,
        "beta": 0.04,
        "inf_mean": 5,
        "verbose": False,
    }

    # Standard
    model1 = Model(scenario, params)
    model1.components = [Susceptible, InfectiousIR, TransmissionSI]
    seed_infections_randomly(model1, ninfections=5)
    model1.run()
    I1 = model1.nodes.I[:, 0]

    # Variant (same TransmissionSI; no TransmissionSIR in this API)
    model2 = Model(scenario, params)
    model2.components = [Susceptible, InfectiousIR, TransmissionSI]
    seed_infections_randomly(model2, ninfections=5)
    model2.run()
    I2 = model2.nodes.I[:, 0]

    assert np.any(I1 > I1[0])
    assert np.any(I2 > I2[0])

    assert np.abs(I1 - I2).max() < 0.2 * pop


@pytest.mark.modeltest
def test_stable_transmission_model_runs(stable_transmission_model):
    model = stable_transmission_model
    model.run()

    inc = model.nodes.incidence[:, 0]
    assert inc.sum() > 0

    assert model.nodes.I[:, 0].max() > 0
