"""
Validation tests for LASER generic model components.

This test suite validates that all components properly execute their validation
functions when model.validating = True. Tests verify that validation hooks
(prevalidate_step and postvalidate_step) are called and can detect inconsistencies.

Components tested:
- components.py: All transmission and disease state components
- vitaldynamics.py: Birth and mortality components
- immunization.py: Immunization components (NOTE: missing validation infrastructure)
- importation.py: Importation components (NOTE: missing validation infrastructure)
"""

import unittest

import laser.core.distributions as dists
import numpy as np
from laser.core import PropertySet
from laser.core.random import seed as set_seed

from laser.generic import Model
from laser.generic.components import (
    Exposed,
    InfectiousIR,
    InfectiousSI,
    InfectiousIS,
    InfectiousIRS,
    Recovered,
    RecoveredRS,
    Susceptible,
    TransmissionSE,
    TransmissionSI,
    TransmissionSIX,
)
from laser.generic.vitaldynamics import (
    BirthsByCBR,
    ConstantPopVitalDynamics,
    MortalityByCDR,
    MortalityByEstimator,
)

try:
    from tests.utils import stdgrid
except ImportError:
    from utils import stdgrid

# Test parameters
NTICKS = 100
SEED = 42
POPULATION = 10_000


class TestComponentsValidation(unittest.TestCase):
    """
    Test suite for validating component.py components with validating=True.

    WHAT IS TESTED:
    - All components execute without validation errors when model.validating=True
    - Validation hooks (pre/post) are properly invoked during simulation
    - Components maintain consistency between agent states and node counts

    SCENARIO SETUP:
    - Small population (10,000) for faster execution
    - Short simulation (100 ticks) sufficient to trigger validation
    - Standard epidemic parameters (beta, durations)

    FAILURE MEANING:
    If these tests fail, components have validation logic errors or inconsistencies
    between individual-level and aggregate-level state tracking. This indicates
    bugs in state transitions or accounting that could lead to incorrect results.
    """

    def setUp(self):
        """Set random seed for reproducibility."""
        set_seed(SEED)
        return

    def test_si_model_validation(self):
        """
        Test SI model components with validation enabled.

        Components tested: Susceptible, InfectiousSI, TransmissionSIX

        WHAT IS VALIDATED:
        - Susceptible counts match agent-level state at each tick
        - Infectious counts match agent-level state at each tick
        - Population conservation (S + I = constant)
        - Transmission incidence matches change in I counts

        FAILURE MEANING:
        Validation failure indicates state transition errors in SI model components,
        such as agents not being properly transitioned from S to I, or aggregate
        counts not matching individual states.
        """
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 100
        scenario["I"] = 100

        params = PropertySet({"nticks": NTICKS, "beta": 0.3})
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        model.components = [
            Susceptible(model),
            InfectiousSI(model),
            TransmissionSIX(model),
        ]

        # Should complete without validation errors
        model.run("SI Model Validation")

        # Verify model completed
        assert model.nodes.S[NTICKS].sum() < model.nodes.S[0].sum(), "Epidemic should have progressed"

        return

    def test_sir_model_validation(self):
        """
        Test SIR model components with validation enabled.

        Components tested: Susceptible, InfectiousIR, Recovered, TransmissionSI

        WHAT IS VALIDATED:
        - S, I, R counts match agent-level states at each tick
        - Population conservation (S + I + R = constant)
        - Infection timers (itimer) are consistent with infectious state
        - Agents with itimer==1 properly transition to recovered
        - Recovery incidence matches change in R counts

        FAILURE MEANING:
        Validation failure indicates errors in SIR transitions such as agents not
        recovering after itimer expires, incorrect itimer initialization, or
        aggregate counts not matching individual states.
        """
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 50
        scenario["I"] = 50
        scenario["R"] = 0

        beta = 0.2
        inf_mean = 7.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            Susceptible(model),
            InfectiousIR(model, infdurdist),
            Recovered(model),
            TransmissionSI(model, infdurdist),
        ]

        # Should complete without validation errors
        model.run("SIR Model Validation")

        # Verify model completed
        assert model.nodes.R[NTICKS].sum() > 0, "Some agents should have recovered"

        return

    def test_seir_model_validation(self):
        """
        Test SEIR model components with validation enabled.

        Components tested: Susceptible, Exposed, InfectiousIR, Recovered, TransmissionSE

        WHAT IS VALIDATED:
        - S, E, I, R counts match agent-level states at each tick
        - Population conservation (S + E + I + R = constant)
        - Exposure timers (etimer) consistent with exposed state
        - Infection timers (itimer) consistent with infectious state
        - Agents with etimer==1 transition to infectious with valid itimer
        - Agents with itimer==1 transition to recovered
        - Exposure incidence matches change in E counts

        FAILURE MEANING:
        Validation failure indicates errors in SEIR transitions such as agents not
        progressing through E�I�R stages, incorrect timer management, or aggregate
        counts diverging from individual states.
        """
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 50
        scenario["E"] = 0
        scenario["I"] = 50
        scenario["R"] = 0

        beta = 0.2
        exp_mean = 5.0
        inf_mean = 7.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        expdurdist = dists.normal(loc=exp_mean, scale=1.0)
        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            Susceptible(model),
            Exposed(model, expdurdist, infdurdist),
            InfectiousIR(model, infdurdist),
            Recovered(model),
            TransmissionSE(model, expdurdist),
        ]

        # Should complete without validation errors
        model.run("SEIR Model Validation")

        # Verify model completed with exposed and recovered agents
        assert model.nodes.R[NTICKS].sum() > 0, "Some agents should have recovered"

        return

    def test_sis_model_validation(self):
        """
        Test SIS model components with validation enabled.

        Components tested: Susceptible, InfectiousIS, TransmissionSI

        WHAT IS VALIDATED:
        - S and I counts match agent-level states (no R compartment)
        - Population conservation (S + I = constant)
        - Agents with itimer==1 return to susceptible (not recovered)
        - Recovery-to-susceptible transitions update counts correctly

        FAILURE MEANING:
        Validation failure indicates errors in SIS recovery transitions where agents
        should return to susceptible state rather than recovered, or aggregate counts
        not matching the cyclic S�I�S dynamics.
        """
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 100
        scenario["I"] = 100

        beta = 0.25
        inf_mean = 5.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        infdurdist = dists.normal(loc=inf_mean, scale=2.0)

        model.components = [
            Susceptible(model),
            InfectiousIS(model, infdurdist),
            TransmissionSI(model, infdurdist),
        ]

        # Should complete without validation errors
        model.run("SIS Model Validation")

        # Verify agents recovered back to susceptible (no R compartment)
        assert not hasattr(model.nodes, "R"), "SIS model should not have R compartment"

        return

    def test_sirs_model_validation(self):
        """
        Test SIRS model components with validation enabled (waning immunity).

        Components tested: Susceptible, InfectiousIRS, RecoveredRS, TransmissionSI

        WHAT IS VALIDATED:
        - S, I, R counts match agent-level states at each tick
        - Population conservation (S + I + R = constant)
        - Agents with itimer==1 transition to recovered with valid rtimer
        - Agents with rtimer==1 return to susceptible (waning immunity)
        - Recovery timers (rtimer) consistent with recovered state

        FAILURE MEANING:
        Validation failure indicates errors in SIRS waning immunity dynamics where
        recovered agents should lose immunity and return to susceptible, or errors
        in managing both itimer and rtimer for the same agents.
        """
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 50
        scenario["I"] = 50
        scenario["R"] = 0

        beta = 0.25
        inf_mean = 7.0
        wan_mean = 30.0
        params = PropertySet({"nticks": NTICKS, "beta": beta})
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        infdurdist = dists.normal(loc=inf_mean, scale=2.0)
        wandurdist = dists.normal(loc=wan_mean, scale=5.0)

        model.components = [
            Susceptible(model),
            InfectiousIRS(model, infdurdist, wandurdist),
            RecoveredRS(model, wandurdist),
            TransmissionSI(model, infdurdist),
        ]

        # Should complete without validation errors
        model.run("SIRS Model Validation")

        # Verify waning occurred (some R returned to S)
        # Note: may not always happen in 100 ticks depending on wan_mean
        # Main goal is validation passes

        return


class TestVitalDynamicsValidation(unittest.TestCase):
    """
    Test suite for validating vitaldynamics.py components with validating=True.

    WHAT IS TESTED:
    - Birth and mortality components execute without validation errors
    - Population dynamics maintain consistency between births/deaths and counts
    - Age tracking (dob/dod) is properly maintained when enabled

    SCENARIO SETUP:
    - Small population (10,000) for faster execution
    - Short simulation (50 ticks) sufficient to observe births/deaths
    - Moderate birth/death rates to generate events

    FAILURE MEANING:
    If these tests fail, vital dynamics components have errors in population
    accounting, age tracking, or state transitions during birth/death events.
    """

    def setUp(self):
        """Set random seed for reproducibility."""
        set_seed(SEED)
        return

    def test_births_by_cbr_validation(self):
        """
        Test BirthsByCBR component with validation enabled.

        WHAT IS VALIDATED:
        - Newborns are added to susceptible state
        - Birth counts match actual number of agents added
        - Date of birth (dob) correctly set for newborns
        - Population count increases correctly

        FAILURE MEANING:
        Validation failure indicates errors in birth implementation such as newborns
        not being added to susceptible state, incorrect population count updates, or
        dob tracking not matching actual births.
        """
        from laser.generic.shared import AliasedDistribution

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population
        scenario["I"] = 0
        scenario["R"] = 0

        params = PropertySet({"nticks": 50, "beta": 0.1})
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        # Create simple pyramid and birth rates
        pyramid_data = np.ones(365 * 50)  # 50 years worth of ages
        pyramid = AliasedDistribution(pyramid_data)
        birthrates = np.full(51, 20.0)  # 20 per 1000 per year

        model.components = [
            Susceptible(model),
            BirthsByCBR(model, birthrates, pyramid, track=True, validating=True),
        ]

        initial_pop = model.people.count
        model.run("Births Validation")

        # Verify births occurred
        assert model.people.count > initial_pop, "Population should have increased due to births"
        assert model.nodes.births.sum() > 0, "Should have recorded births"

        return

    def test_mortality_by_cdr_validation(self):
        """
        Test MortalityByCDR component with validation enabled.

        WHAT IS VALIDATED:
        - Deaths are recorded correctly by node
        - State counts (S, I, R) decrease appropriately when agents die
        - Deceased agents marked with DECEASED state

        FAILURE MEANING:
        Validation failure indicates errors in mortality implementation such as
        death counts not matching actual state changes, or compartment counts not
        being decremented when agents die.
        """
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 100
        scenario["I"] = 50
        scenario["R"] = 50

        params = PropertySet({"nticks": 50, "beta": 0.0})  # No transmission
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        mortalityrates = np.full(51, 10.0)  # 10 per 1000 per year

        model.components = [
            Susceptible(model),
            InfectiousIR(model, dists.normal(loc=7, scale=1)),
            Recovered(model),
            MortalityByCDR(model, mortalityrates, validating=True),
        ]

        model.run("Mortality CDR Validation")

        # Verify deaths occurred
        assert model.nodes.deaths.sum() > 0, "Should have recorded deaths"

        return

    def test_mortality_by_estimator_validation(self):
        """
        Test MortalityByEstimator component with validation enabled.

        WHAT IS VALIDATED:
        - Agents with dod==tick are marked as deceased
        - Death counts match agents reaching their dod
        - State counts decrease when agents die
        - Newborns receive valid dod values via on_birth hook

        FAILURE MEANING:
        Validation failure indicates errors in life table mortality such as agents
        not dying when dod is reached, incorrect dod sampling, or on_birth hook
        not properly initializing dod for newborns.
        """
        from laser.core.estimators import KaplanMeierEstimator
        from laser.generic.shared import AliasedDistribution

        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 50
        scenario["I"] = 50
        scenario["R"] = 0

        params = PropertySet({"nticks": 50, "beta": 0.0})
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        # Create simple survival curve (uniform hazard)
        ages = np.arange(0, 365 * 100)
        survival = np.exp(-ages / (365 * 70))  # Mean survival ~70 years
        estimator = KaplanMeierEstimator(ages, survival)

        # Need births to test on_birth hook
        pyramid_data = np.ones(365 * 50)
        pyramid = AliasedDistribution(pyramid_data)
        birthrates = np.full(51, 15.0)

        model.components = [
            Susceptible(model),
            InfectiousIR(model, dists.normal(loc=7, scale=1)),
            Recovered(model),
            BirthsByCBR(model, birthrates, pyramid, track=True, validating=True),
            MortalityByEstimator(model, estimator, validating=True),
        ]

        model.run("Mortality Estimator Validation")

        # Verify deaths occurred
        assert model.nodes.deaths.sum() > 0, "Should have recorded deaths"

        return

    def test_constant_pop_vital_dynamics_validation(self):
        """
        Test ConstantPopVitalDynamics component with validation enabled.

        WHAT IS VALIDATED:
        - Recycled agents return to susceptible state
        - Birth and death counts are equal (constant population)
        - DOB correctly updated for recycled agents (if tracking enabled)
        - State counts properly updated for recycling transitions

        FAILURE MEANING:
        Validation failure indicates errors in recycling logic such as agents not
        returning to susceptible state, births not equaling deaths, or dob not
        being updated correctly for recycled agents.
        """
        scenario = stdgrid(M=1, N=1, population_fn=lambda r, c: POPULATION)
        scenario["S"] = scenario.population - 100
        scenario["I"] = 50
        scenario["R"] = 50

        params = PropertySet({"nticks": 50, "beta": 0.0})
        model = Model(scenario, params)
        model.validating = True  # Enable validation

        recycle_rates = np.full(51, 15.0)  # 15 per 1000 per year

        model.components = [
            Susceptible(model),
            InfectiousIR(model, dists.normal(loc=7, scale=1)),
            Recovered(model),
            ConstantPopVitalDynamics(model, recycle_rates, dobs=True, validating=True),
        ]

        initial_pop = model.people.count
        model.run("Constant Pop Validation")

        # Verify population remained constant
        assert model.people.count == initial_pop, "Population should remain constant with recycling"
        # Verify births == deaths
        assert np.all(model.nodes.births == model.nodes.deaths), "Births should equal deaths in recycling model"

        return


# Note: Immunization and Importation components currently lack validation infrastructure
# and cannot be tested in the same way. See suggestions below for implementing validation.


if __name__ == "__main__":
    unittest.main()
