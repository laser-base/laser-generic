from typing import Callable
from typing import Union

import laser.core.distributions as dists
import matplotlib.pyplot as plt
import numba as nb
import numpy as np

from laser.generic.shared import State
from laser.generic.utils import ValuesMap
from laser.generic.utils import validate


@nb.njit(
    nogil=True,
    parallel=True,
    cache=True,
)
def nb_timer_update(states, test_state, timers, new_state, transitioned, node_ids):
    for i in nb.prange(len(states)):
        if states[i] == test_state:
            timers[i] -= 1
            if timers[i] == 0:
                states[i] = new_state
                transitioned[nb.get_thread_id(), node_ids[i]] += 1

    return


@nb.njit(
    nogil=True,
    parallel=True,
)
def nb_timer_update_timer_set(
    states, test_state, oldtimers, new_state, newtimers, transitioned, node_ids, duration_dist, duration_min, tick
):
    for i in nb.prange(len(states)):
        if states[i] == test_state:
            oldtimers[i] -= 1
            if oldtimers[i] == 0:
                states[i] = new_state
                nid = node_ids[i]
                newtimers[i] = np.maximum(np.round(duration_dist(tick, nid)), duration_min)  # Set the new timer
                transitioned[nb.get_thread_id(), nid] += 1

    return


class Susceptible:
    """
    Susceptible Component for Patch-Based Agent-Based Models (S, SI, SIS, SIR, SEIR, etc.)

    This component initializes and tracks the count of susceptible individuals (`S`) in
    a spatially structured agent-based model. It is compatible with all standard LASER
    disease progression models that include a "susceptible" state.

    Responsibilities:
    - Initializes agent-level properties:
        • `nodeid`: Patch ID of each agent (uint16)
        • `state`: Infection state (int8), defaulting to `State.SUSCEPTIBLE`
    - Initializes node-level property:
        • `S[t, i]`: Susceptible count in node `i` at time `t`
    - At each timestep, propagates the susceptible count forward (`S[t+1] = S[t]`),
      unless modified by other components (e.g., exposure, births).
    - Validates consistency between patch-level susceptible counts and agent-level state.

    Usage:
    Add this component early in the component list for any model with SUSCEPTIBLE agents,
    typically before transmission or exposure components. Compatible with:
        - `SIR.Transmission`
        - `SIR.Exposure`
        - `SIR.Infectious`
        - `SIR.Recovered`
        - Custom SEIRS extensions

    Requires:
    - `model.people`: A LaserFrame for all agents
    - `model.nodes`: Patch-level state
    - `model.scenario`: Input DataFrame with `population` and optionally `S` columns
    - `model.params.nticks`: Number of simulation ticks

    Validation:
    - Ensures consistency of susceptible counts before and after each step
    - Prevents unintentional state drift by validating against agent `state` values

    Output:
    - `model.nodes.S`: A `(nticks+1, num_nodes)` array of susceptible counts
    - Optional plotting via `plot()` for visual inspection of per-node and total `S`

    Step Behavior:
        For tick t:
            S[t+1] = S[t]     # Unless explicitly modified by other components

    This component does not alter agent states directly but serves as a synchronized
    counter and validator of susceptible individuals.

    Example:
        model.components = [
            SIR.Susceptible(model),
            SIR.Transmission(model, ...),
            SIR.Exposure(model),
            SIR.Infectious(model, ...),
            SIR.Recovered(model),
        ]
    """

    def __init__(self, model, validating: bool = False):
        self.model = model
        self.validating = validating

        self.model.people.add_scalar_property("nodeid", dtype=np.uint16)
        self.model.people.add_scalar_property("state", dtype=np.int8, default=State.SUSCEPTIBLE.value)
        self.model.nodes.add_vector_property("S", model.params.nticks + 1, dtype=np.int32)

        self.model.people.nodeid[:] = np.repeat(np.arange(self.model.nodes.count, dtype=np.uint16), self.model.scenario.population)
        self.model.nodes.S[0] = self.model.scenario.S

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_unchanged(self.model.nodes.S[tick], self.model.nodes.S[tick + 1], "Susceptible counts")

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        return

    def plot(self):
        _fig, ax1 = plt.subplots()
        for node in range(self.model.nodes.count):
            ax1.plot(self.model.nodes.S[:, node], label=f"Node {node}")
        ax1.set_xlabel("Tick")
        ax1.set_ylabel("Susceptible (by Node)")
        ax1.set_title("Susceptible over Time by Node")
        ax1.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(np.sum(self.model.nodes.S, axis=1), color="black", linestyle="--", label="Total Susceptible")
        ax2.set_ylabel("Total Susceptible")
        ax2.legend(loc="upper right")

        plt.show()

        return


class Exposed:
    """
    Exposed Component for SEIR/SEIRS Models with Explicit Incubation Period

    This component handles the incubation phase in models where agents must transition from
    an 'exposed' (E) state to 'infectious' (I) after a delay. It supports custom incubation
    and infectious duration distributions and handles both initialization and per-tick dynamics.

    Agents transition from Exposed to Infectious when their incubation timer (etimer) expires.
    Tracks number of agents becoming infectious each tick in `model.nodes.newly_infectious`.

    Responsibilities:
    - Initializes exposed individuals at time 0 (if provided in the scenario)
    - Assigns and tracks per-agent incubation timers (`etimer`)
    - Transitions agents from `EXPOSED` to `INFECTIOUS` when `etimer == 0`
    - Assigns new infection timers (`itimer`) upon becoming infectious
    - Updates patch-level EXPOSED (`E`) and INFECTIOUS case counts
    - Provides validation hooks for state and timer consistency

    Required Inputs:
    - `model.scenario.E`: initial count of exposed individuals per node (optional)
    - `expdurdist`: callable returning sampled incubation durations
    - `infdurdist`: callable returning sampled infectious durations
    - `expdurmin`: minimum incubation period (default 1 day)
    - `infdurmin`: minimum infectious period (default 1 day)

    Outputs:
    - `model.people.etimer`: agent-level incubation timer
    - `model.nodes.E[t, i]`: number of exposed individuals at time `t` in node `i`
    - `model.nodes.newly_infectious[t, i]`: number of newly infectious cases per node per day

    Validation:
    - Ensures consistency between individual states and `etimer` values
    - Ensures that agents becoming infectious have valid `itimer` values assigned
    - Prevents agents with expired `etimer` from remaining in EXPOSED state

    Step Behavior:
        For each agent:
            - Decrease `etimer`
            - If `etimer == 0`, change state to `INFECTIOUS` and assign `itimer`
            - Update `model.nodes.E` and `model.nodes.I` counts accordingly

    Plotting:
    The `plot()` method provides a time series of exposed individuals per node and total across all nodes.

    Example:
        model.components = [
            SIR.Susceptible(model),
            Exposed(model, expdurdist, infdurdist),
            SIR.Infectious(model, infdurdist),
            ...
        ]
    """

    def __init__(self, model, expdurdist, infdurdist, expdurmin=1, infdurmin=1, validating: bool = False):
        self.model = model
        self.validating = validating

        self.model.people.add_scalar_property("etimer", dtype=np.uint16)
        self.model.nodes.add_vector_property("E", model.params.nticks + 1, dtype=np.int32)
        self.model.nodes.add_vector_property("newly_infectious", model.params.nticks + 1, dtype=np.int32)

        self.model.nodes.E[0] = self.model.scenario.E

        self.expdurdist = expdurdist
        self.infdurdist = infdurdist
        self.expdurmin = expdurmin
        self.infdurmin = infdurmin

        # convenience
        nodeids = self.model.people.nodeid
        states = self.model.people.state

        for node in range(self.model.nodes.count):
            nseeds = self.model.scenario.E[node]
            if nseeds > 0:
                i_susceptible = np.nonzero((nodeids == node) & (states == State.SUSCEPTIBLE.value))[0]
                assert nseeds <= len(i_susceptible), (
                    f"Node {node} has more initial exposed ({nseeds}) than available susceptible ({len(i_susceptible)})"
                )
                i_exposed = np.random.choice(i_susceptible, size=nseeds, replace=False)
                self.model.people.state[i_exposed] = State.EXPOSED.value
                samples = dists.sample_floats(self.expdurdist, np.zeros(nseeds, dtype=np.float32))
                samples = np.round(samples)
                samples = np.maximum(samples, self.expdurmin).astype(self.model.people.etimer.dtype)
                self.model.people.etimer[i_exposed] = samples
                assert np.all(self.model.people.etimer[i_exposed] > 0), (
                    f"Exposed individuals should have etimer > 0 ({self.model.people.etimer[i_exposed].min()=})"
                )

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.E[tick + 1], self.model.people, State.EXPOSED, "Exposed")
        _check_timer_active(self.model.people.state, State.EXPOSED.value, self.model.people.etimer, "Exposed", "etimer")
        _check_state_timer_consistency(self.model.people.state, State.EXPOSED.value, self.model.people.etimer, "Exposed", "etimer")

        alive = self.model.people.state != State.DECEASED.value
        self.etimers_one = (alive) & (self.model.people.etimer == 1)

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.E[tick + 1], self.model.people, State.EXPOSED, "Exposed")
        _check_timer_active(self.model.people.state, State.EXPOSED.value, self.model.people.etimer, "Exposed", "etimer")
        _check_state_timer_consistency(self.model.people.state, State.EXPOSED.value, self.model.people.etimer, "Exposed", "etimer")

        assert np.all(self.model.people.state[self.etimers_one] == State.INFECTIOUS.value), (
            "Individuals with etimer == 1 before should now be infectious."
        )
        assert np.all(self.model.people.etimer[self.etimers_one] == 0), "Individuals with etimer == 1 before should now have etimer == 0."
        assert np.all(self.model.people.itimer[self.etimers_one] > 0), "Individuals with etimer == 1 before should now have itimer > 0."

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        newly_infectious_by_node = np.zeros((nb.get_num_threads(), self.model.nodes.count), dtype=np.int32)
        nb_timer_update_timer_set(
            self.model.people.state,
            State.EXPOSED.value,
            self.model.people.etimer,
            State.INFECTIOUS.value,
            self.model.people.itimer,
            newly_infectious_by_node,
            self.model.people.nodeid,
            self.infdurdist,
            self.infdurmin,
            tick,
        )
        newly_infectious_by_node = newly_infectious_by_node.sum(axis=0).astype(self.model.nodes.S.dtype)  # Sum over threads

        # state(t+1) = state(t) + ∆state(t)
        self.model.nodes.E[tick + 1] -= newly_infectious_by_node
        self.model.nodes.I[tick + 1] += newly_infectious_by_node
        # Record today's ∆
        self.model.nodes.newly_infectious[tick] = newly_infectious_by_node

        return

    def plot(self):
        _fig, ax1 = plt.subplots()
        for node in range(self.model.nodes.count):
            ax1.plot(self.model.nodes.E[:, node], label=f"Node {node}")
        ax1.set_xlabel("Tick")
        ax1.set_ylabel("Exposed (by Node)")
        ax1.set_title("Exposed over Time by Node")
        ax1.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(np.sum(self.model.nodes.E, axis=1), color="black", linestyle="--", label="Total Exposed")
        ax2.set_ylabel("Total Exposed")
        ax2.legend(loc="upper right")
        plt.show()

        return


class InfectiousSI:
    """
    Infectious Component for SI Models (No Recovery)

    This component manages the infectious state in SI-style epidemic models where agents
    remain infectious indefinitely. It is appropriate for use in models without a recovered
    or removed state (i.e., no `R` compartment).

    Responsibilities:
    - Initializes agents as infectious based on `model.scenario.I`
    - Tracks the number of infectious individuals (`I`) in each patch over time
    - Maintains per-tick, per-node counts in `model.nodes.I`
    - Validates consistency between agent states and patch-level totals

    Required Inputs:
    - `model.scenario.I`: array of initial infected counts per patch
    - `model.people.state`: infection state per agent
    - `model.people.nodeid`: patch assignment per agent
    - `model.params.nticks`: number of timesteps to simulate

    Outputs:
    - `model.nodes.I[t, i]`: number of infectious individuals in node `i` at time `t`

    Step Behavior:
        For each timestep `t`, this component copies:
            I[t+1] = I[t]
        (No recovery or removal; new infections may be added externally.)

    Validation:
    - Ensures that patch-level infectious counts (`model.nodes.I`) match the agent-level state
    - Asserts that the sum of `S` and `I` matches total population at initialization
    - Validates that infected counts do not change unexpectedly (unless altered by another component)

    Plotting:
    The `plot()` method shows the number of infectious agents per patch and in total across time.

    Example:
        model.components = [
            SIR.Susceptible(model),
            InfectiousSI(model),
            SIR.Transmission(model, ...),
        ]
    """

    def __init__(self, model, validating: bool = False):
        self.model = model
        self.validating = validating

        self.model.nodes.add_vector_property("I", model.params.nticks + 1, dtype=np.int32)

        # convenience
        nodeids = self.model.people.nodeid
        populations = self.model.scenario.population.values

        for node in range(self.model.nodes.count):
            citizens = np.nonzero(nodeids == node)[0]
            assert len(citizens) == populations[node], f"Found {len(citizens)} citizens in node {node} but expected {populations[node]}"
            nseeds = self.model.scenario.I[node]
            assert nseeds <= len(citizens), f"Node {node} has more initial infected ({nseeds}) than population ({len(citizens)})"
            if nseeds > 0:
                indices = np.random.choice(citizens, size=nseeds, replace=False)
                self.model.people.state[indices] = State.INFECTIOUS.value

        self.model.nodes.I[0] = self.model.scenario.I
        assert np.all(self.model.nodes.S[0] + self.model.nodes.I[0] == self.model.scenario.population), (
            "Initial S + I does not equal population"
        )

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")
        _check_unchanged(self.model.nodes.I[tick], self.model.nodes.I[tick + 1], "Infected counts")

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        """Step function for the Infected component.

        Args:
            tick (int): The current tick of the simulation.
        """

        return

    def plot(self):
        _fig, ax1 = plt.subplots()
        for node in range(self.model.nodes.count):
            ax1.plot(self.model.nodes.I[:, node], label=f"Node {node}")
        ax1.set_xlabel("Tick")
        ax1.set_ylabel("Infected (by Node)")
        ax1.set_title("Infected over Time by Node")
        ax1.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(np.sum(self.model.nodes.I, axis=1), color="black", linestyle="--", label="Total Infected")
        ax2.set_ylabel("Total Infected")
        ax2.legend(loc="upper right")

        plt.show()

        return


class InfectiousIS:
    """
    Infectious Component for SIS Models (Infection + Recovery to Susceptible)

    This component handles the infectious state in SIS-style models, where agents
    recover from infection and immediately return to the susceptible pool. It supports
    per-agent infection durations and manages patch-level infectious counts over time.

    Agents transition from Infectious back to Susceptible after the infectious period (itimer).
    Tracks number of agents recovering each tick in `model.nodes.newly_recovered`.

    Responsibilities:
    - Initializes infected agents and their infection timers (`itimer`)
    - Decrements `itimer` daily for infectious agents
    - Automatically transitions agents from `INFECTIOUS` to `SUSCEPTIBLE` when `itimer == 0`
    - Tracks per-day recoveries at the node level in `model.nodes.recovered`
    - Maintains node-level `I` and `S` counts with full timestep resolution

    Required Inputs:
    - `model.scenario.I`: initial number of infectious agents per patch
    - `infdurdist`: a callable function which samples the infectious duration distribution
    - `infdurmin`: the minimum infection period (default = 1 time step)

    Outputs:
    - `model.people.itimer`: per-agent infection countdown timer
    - `model.nodes.I[t, i]`: number of infectious individuals at tick `t` in node `i`
    - `model.nodes.recovered[t, i]`: number of recoveries at tick `t` in node `i`

    Step Behavior:
        At each tick:
        - Infectious agents decrement their `itimer`
        - Agents with `itimer == 0` are transitioned back to susceptible
        - `model.nodes.I` is updated accordingly
        - Recovered counts are recorded in `model.nodes.recovered`

    Validation:
    - Ensures consistency between agent `state` and infection timer (`itimer`)
    - Validates `I` census against agent-level state before and after each tick

    Plotting:
    The `plot()` method displays both per-node and total infectious counts over time.

    Example:
        model.components = [
            SIR.Susceptible(model),
            InfectiousIS(model, infdurdist),
            SIR.Transmission(model, ...),
        ]
    """

    def __init__(self, model, infdurdist, infdurmin=1, validating: bool = False):
        self.model = model
        self.validating = validating

        self.model.people.add_scalar_property("itimer", dtype=np.uint16)
        self.model.nodes.add_vector_property("I", model.params.nticks + 1, dtype=np.int32)
        self.model.nodes.add_vector_property("newly_recovered", model.params.nticks + 1, dtype=np.int32)

        self.infdurdist = infdurdist
        self.infdurmin = infdurmin

        # convenience
        nodeids = self.model.people.nodeid
        populations = self.model.scenario.population.values

        for node in range(self.model.nodes.count):
            nseeds = self.model.scenario.I[node]
            if nseeds > 0:
                i_citizens = np.nonzero(nodeids == node)[0]
                assert len(i_citizens) == populations[node], (
                    f"Found {len(i_citizens)} citizens in node {node} but expected {populations[node]}"
                )
                assert nseeds <= len(i_citizens), f"Node {node} has more initial infectious ({nseeds}) than population ({len(i_citizens)})"
                i_infectious = np.random.choice(i_citizens, size=nseeds, replace=False)
                self.model.people.state[i_infectious] = State.INFECTIOUS.value
                samples = dists.sample_floats(self.infdurdist, np.zeros(nseeds, dtype=np.int32))
                samples = np.round(samples)
                samples = np.maximum(samples, self.infdurmin).astype(self.model.people.itimer.dtype)
                self.model.people.itimer[i_infectious] = samples
                assert np.all(self.model.people.itimer[i_infectious] > 0), (
                    f"Infected individuals should have itimer > 0 ({self.model.people.itimer[i_infectious].min()=})"
                )

        self.model.nodes.I[0] = self.model.scenario.I
        assert np.all(self.model.nodes.S[0] + self.model.nodes.I[0] == self.model.scenario.population), (
            "Initial S + I does not equal population"
        )

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")
        _check_timer_active(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")
        _check_state_timer_consistency(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")
        _check_timer_active(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")
        _check_state_timer_consistency(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        """Step function for the Infected component.

        Args:
            tick (int): The current tick of the simulation.
        """

        newly_recovered_by_node = np.zeros((nb.get_num_threads(), self.model.nodes.count), dtype=np.int32)
        nb_timer_update(
            self.model.people.state,
            State.INFECTIOUS.value,
            self.model.people.itimer,
            State.SUSCEPTIBLE.value,
            newly_recovered_by_node,
            self.model.people.nodeid,
        )
        newly_recovered_by_node = newly_recovered_by_node.sum(axis=0).astype(self.model.nodes.S.dtype)  # Sum over threads

        # state(t+1) = state(t) + ∆state(t)
        self.model.nodes.S[tick + 1] += newly_recovered_by_node
        self.model.nodes.I[tick + 1] -= newly_recovered_by_node
        # Record today's ∆
        self.model.nodes.newly_recovered[tick] = newly_recovered_by_node

        return

    def plot(self):
        _fig, ax1 = plt.subplots()
        for node in range(self.model.nodes.count):
            ax1.plot(self.model.nodes.I[:, node], label=f"Node {node}")
        ax1.set_xlabel("Tick")
        ax1.set_ylabel("Infected (by Node)")
        ax1.set_title("Infected over Time by Node")
        ax1.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(np.sum(self.model.nodes.I, axis=1), color="black", linestyle="--", label="Total Infected")
        ax2.set_ylabel("Total Infected")
        ax2.legend(loc="upper right")
        plt.show()

        return


class InfectiousIR:
    """
    Infectious Component for SIR/SEIR Models (With Recovery to Immune)

    This component manages agents in the infectious state for models where infected individuals
    recover permanently (i.e., transition to a `RECOVERED` state without waning). It supports
    agent-level infection durations and patch-level tracking of recoveries over time.

    Infectious component for an SIR/SEIR model - includes infectious duration, no waning immunity in newly_recovered state.

    Agents transition from Infectious to Recovered after the infectious period (itimer).
    Tracks number of agents recovering each tick in `model.nodes.newly_recovered`.

    Responsibilities:
    - Initializes infected agents and their infection timers (`itimer`) based on scenario input
    - Decrements `itimer` daily for infectious agents
    - Transitions agents from `INFECTIOUS` to `RECOVERED` when `itimer == 0`
    - Updates patch-level state variables:
        • `I[t, i]`: infectious count at tick `t` in node `i`
        • `R[t, i]`: recovered count
        • `recovered[t, i]`: number of recoveries during tick `t`

    Required Inputs:
    - `model.scenario.I`: number of initially infected individuals per patch
    - `infdurdist`: function returning infection durations
    - `infdurmin`: minimum infectious period (default = 1 day)

    Outputs:
    - `model.people.itimer`: countdown timers per agent
    - `model.nodes.I[t]`, `.R[t]`: infectious and recovered counts per patch
    - `model.nodes.recovered[t]`: daily recoveries per patch

    Step Behavior:
    - Infectious agents decrement `itimer`
    - When `itimer == 0`, agent state is set to `RECOVERED`
    - Patch-level `I` and `R` are updated; `recovered` logs today's transitions

    Validation:
    - Ensures internal consistency between agent state and timer
    - Confirms agents with `itimer == 1` recover exactly one day later
    - Validates population conservation (`S + I + R = N`)

    Plotting:
    The `plot()` method shows per-node and total infectious counts across time.

    Example:
        model.components = [
            SIR.Susceptible(model),
            InfectiousIR(model, infdurdist),
            SIR.Recovered(model),
            ...
        ]
    """

    def __init__(self, model, infdurdist, infdurmin=1, validating: bool = False):
        self.model = model
        self.validating = validating

        self.model.people.add_scalar_property("itimer", dtype=np.uint16)
        self.model.nodes.add_vector_property("I", model.params.nticks + 1, dtype=np.int32)
        self.model.nodes.add_vector_property("newly_recovered", model.params.nticks + 1, dtype=np.int32)

        self.model.nodes.I[0] = self.model.scenario.I

        self.infdurdist = infdurdist
        self.infdurmin = infdurmin

        # convenience
        nodeids = self.model.people.nodeid
        states = self.model.people.state

        for node in range(self.model.nodes.count):
            nseeds = self.model.scenario.I[node]
            if nseeds > 0:
                i_susceptible = np.nonzero((nodeids == node) & (states == State.SUSCEPTIBLE.value))[0]
                assert nseeds <= len(i_susceptible), (
                    f"Node {node} has more initial infectious ({nseeds}) than available susceptible ({len(i_susceptible)})"
                )
                i_infectious = np.random.choice(i_susceptible, size=nseeds, replace=False)
                self.model.people.state[i_infectious] = State.INFECTIOUS.value
                samples = dists.sample_floats(self.infdurdist, np.zeros(nseeds, dtype=np.float32))
                samples = np.round(samples)
                samples = np.maximum(samples, self.infdurmin).astype(self.model.people.itimer.dtype)
                self.model.people.itimer[i_infectious] = samples
                assert np.all(self.model.people.itimer[i_infectious] > 0), (
                    f"Infected individuals should have itimer > 0 ({self.model.people.itimer[i_infectious].min()=})"
                )

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")
        _check_timer_active(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")
        _check_state_timer_consistency(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")

        alive = self.model.people.state != State.DECEASED.value
        self.itimers_one = (alive) & (self.model.people.itimer == 1)

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")
        _check_timer_active(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")
        _check_state_timer_consistency(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")

        assert np.all(self.model.people.state[self.itimers_one] == State.RECOVERED.value), (
            "Individuals with itimer == 1 before should now be recovered."
        )
        assert np.all(self.model.people.itimer[self.itimers_one] == 0), "Individuals with itimer == 1 before should now have itimer == 0."

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        """Step function for the Infected component.

        Args:
            tick (int): The current tick of the simulation.
        """

        newly_recovered_by_node = np.zeros((nb.get_num_threads(), self.model.nodes.count), dtype=np.int32)
        nb_timer_update(
            self.model.people.state,
            State.INFECTIOUS.value,
            self.model.people.itimer,
            State.RECOVERED.value,
            newly_recovered_by_node,
            self.model.people.nodeid,
        )
        newly_recovered_by_node = newly_recovered_by_node.sum(axis=0).astype(self.model.nodes.S.dtype)  # Sum over threads

        # state(t+1) = state(t) + ∆state(t)
        self.model.nodes.I[tick + 1] -= newly_recovered_by_node
        self.model.nodes.R[tick + 1] += newly_recovered_by_node
        # Record today's ∆
        self.model.nodes.newly_recovered[tick] = newly_recovered_by_node

        return

    def plot(self):
        _fig, ax1 = plt.subplots()
        for node in range(self.model.nodes.count):
            ax1.plot(self.model.nodes.I[:, node], label=f"Node {node}")
        ax1.set_xlabel("Tick")
        ax1.set_ylabel("Infected (by Node)")
        ax1.set_title("Infected over Time by Node")
        ax1.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(np.sum(self.model.nodes.I, axis=1), color="black", linestyle="--", label="Total Infected")
        ax2.set_ylabel("Total Infected")
        ax2.legend(loc="upper right")
        plt.show()

        return


class InfectiousIRS:
    """
    Infectious Component for SIRS/SEIRS Models (Recovery with Waning Immunity)

    This component manages infectious individuals in models where recovery confers
    temporary immunity, after which agents become susceptible again (SIRS/SEIRS).

    Agents transition from Infectious to Recovered after the infectious period (itimer).
    Set the waning immunity timer (rtimer) upon recovery.
    Tracks number of agents recovering each tick in `model.nodes.newly_recovered`.

    Responsibilities:
    - Initializes infectious agents from `model.scenario.I`
    - Assigns and tracks infectious timers (`itimer`) per agent
    - Transitions agents from `INFECTIOUS` to `RECOVERED` when `itimer == 0`
    - Assigns a waning immunity timer (`rtimer`) upon recovery
    - Updates patch-level state:
        • `I[t, i]`: current infectious count
        • `R[t, i]`: current recovered count
        • `recovered[t, i]`: number of agents recovering on tick `t`

    Required Inputs:
    - `model.scenario.I`: number of initially infected agents per node
    - `infdurdist`: function that samples the infectious duration distribution
    - `wandurdist`: function that samples the waning immunity duration distribution
    - `infdurmin`: minimum infectious period (default = 1 day)
    - `wandurmin`: minimum duration of immunity (default = 1 day)

    Outputs:
    - `model.people.itimer`: days remaining in the infectious state
    - `model.people.rtimer`: days remaining in the recovered state
    - `model.nodes.I`, `model.nodes.R`: counts per node per tick
    - `model.nodes.recovered[t]`: number of recoveries recorded on tick `t`

    Step Behavior:
    - Infectious agents decrement their `itimer`
    - When `itimer == 0`, agents become recovered and receive an `rtimer`
    - Patch-level totals are updated
    - Downstream components (e.g., `Recovered`) handle `rtimer` countdown and eventual return to `SUSCEPTIBLE`

    Validation:
    - Ensures timer consistency and population accounting
    - Confirms correct infectious-to-recovered transitions
    - Can be chained with recovery and waning components for full SIRS/SEIRS loops

    Plotting:
    Two plots are provided:
    1. Infected counts per node
    2. Total infected and recovered counts across time

    Example:
        model.components = [
            SIR.Susceptible(model),
            InfectiousIRS(model, infdurdist, wandurdist),
            Exposed(model, ...),
            Recovered(model),
        ]
    """

    def __init__(self, model, infdurdist, wandurdist, infdurmin=1, wandurmin=1, validating: bool = False):
        self.model = model
        self.validating = validating

        self.model.people.add_scalar_property("itimer", dtype=np.uint16)
        if not hasattr(self.model.people, "rtimer"):
            self.model.people.add_scalar_property("rtimer", dtype=np.uint16)
        self.model.nodes.add_vector_property("I", model.params.nticks + 1, dtype=np.int32)
        self.model.nodes.add_vector_property("newly_recovered", model.params.nticks + 1, dtype=np.int32)

        self.model.nodes.I[0] = self.model.scenario.I

        self.infdurdist = infdurdist
        self.wandurdist = wandurdist
        self.infdurmin = infdurmin
        self.wandurmin = wandurmin

        # convenience
        nodeids = self.model.people.nodeid
        states = self.model.people.state

        for node in range(self.model.nodes.count):
            nseeds = self.model.scenario.I[node]
            if nseeds > 0:
                i_susceptible = np.nonzero((nodeids == node) & (states == State.SUSCEPTIBLE.value))[0]
                assert nseeds <= len(i_susceptible), (
                    f"Node {node} has more initial infectious ({nseeds}) than available susceptible ({len(i_susceptible)})"
                )
                i_infectious = np.random.choice(i_susceptible, size=nseeds, replace=False)
                self.model.people.state[i_infectious] = State.INFECTIOUS.value
                samples = dists.sample_floats(self.infdurdist, np.zeros(nseeds, dtype=np.float32))
                samples = np.round(samples)
                samples = np.maximum(self.infdurmin, samples).astype(self.model.people.itimer.dtype)
                self.model.people.itimer[i_infectious] = samples
                assert np.all(self.model.people.itimer[i_infectious] > 0), (
                    f"Infected individuals should have itimer > 0 ({self.model.people.itimer[i_infectious].min()=})"
                )

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")
        _check_timer_active(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")
        _check_state_timer_consistency(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")
        _check_timer_active(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")
        _check_state_timer_consistency(self.model.people.state, State.INFECTIOUS.value, self.model.people.itimer, "Infectious", "itimer")

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        """Step function for the Infected component.

        Args:
            tick (int): The current tick of the simulation.
        """

        newly_recovered_by_node = np.zeros((nb.get_num_threads(), self.model.nodes.count), dtype=np.int32)
        nb_timer_update_timer_set(
            self.model.people.state,
            State.INFECTIOUS.value,
            self.model.people.itimer,
            State.RECOVERED.value,
            self.model.people.rtimer,
            newly_recovered_by_node,
            self.model.people.nodeid,
            self.wandurdist,
            self.wandurmin,
            tick,
        )
        newly_recovered_by_node = newly_recovered_by_node.sum(axis=0).astype(self.model.nodes.S.dtype)  # Sum over threads

        # state(t+1) = state(t) + ∆state(t)
        self.model.nodes.I[tick + 1] -= newly_recovered_by_node
        self.model.nodes.R[tick + 1] += newly_recovered_by_node
        # Record today's ∆
        self.model.nodes.newly_recovered[tick] = newly_recovered_by_node

        return

    def plot(self):
        # First plot: Infected over Time by Node
        _fig, ax1 = plt.subplots()
        for node in range(self.model.nodes.count):
            ax1.plot(self.model.nodes.I[:, node], label=f"Node {node}")
        ax1.set_xlabel("Tick")
        ax1.set_ylabel("Infected (by Node)")
        ax1.set_title("Infected over Time by Node")
        ax1.legend(loc="upper left")
        plt.show()

        # Second plot: Total Infected and Total Recovered over Time
        _fig, ax2 = plt.subplots()
        total_infected = np.sum(self.model.nodes.I, axis=1)
        total_recovered = np.sum(self.model.nodes.newly_recovered, axis=1)
        ax2.plot(total_infected, color="black", linestyle="--", label="Total Infected")
        ax2.plot(total_recovered, color="green", linestyle="-.", label="Total Recovered")
        ax2.set_xlabel("Tick")
        ax2.set_ylabel("Count")
        ax2.set_title("Total Infected and Total Recovered Over Time")
        ax2.legend(loc="upper right")
        plt.show()

        return


class Recovered:
    """
    Recovered Component for SIR/SEIR Models (Permanent Immunity)

    This component manages agents in the recovered state in models where immunity does
    not wane (i.e., once recovered, agents stay recovered permanently). It tracks the
    number of recovered individuals over time at the patch level, but performs no active
    transitions itself — recovery transitions must be handled by upstream components.

    Responsibilities:
    - Initializes agents as recovered if specified in `model.scenario.R`
    - Tracks per-patch recovered counts over time in `model.nodes.R`
    - Verifies consistency between agent state and aggregate recovered counts
    - Propagates recovered totals forward unchanged (unless modified by other components)

    Required Inputs:
    - `model.scenario.R`: number of initially recovered individuals per node

    Outputs:
    - `model.nodes.R[t, i]`: number of recovered individuals at tick `t` in node `i`

    Step Behavior:
    - At each tick, carries forward:
        R[t+1] = R[t]
    - This component does not change any agent's state or internal timers

    🧪 Validation:
    - Ensures per-agent state matches aggregate `R` counts before and after each step
    - Detects accidental changes to recovered counts not explained by upstream logic

    Plotting:
    The `plot()` method shows per-node and total recovered counts over time.

    Example:
        model.components = [
            SIR.Susceptible(model),
            InfectiousIR(model, infdurdist),
            Recovered(model),  # passive tracker, assumes recovery handled upstream
        ]
    """

    def __init__(self, model, validating: bool = False):
        self.model = model
        self.validating = validating

        self.model.nodes.add_vector_property("R", model.params.nticks + 1, dtype=np.int32)

        self.model.nodes.R[0] = self.model.scenario.R

        # convenience
        nodeids = self.model.people.nodeid
        states = self.model.people.state

        for node in range(self.model.nodes.count):
            nseeds = self.model.scenario.R[node]
            if nseeds > 0:
                i_susceptible = np.nonzero((nodeids == node) & (states == State.SUSCEPTIBLE.value))[0]
                assert nseeds <= len(i_susceptible), (
                    f"Node {node} has more initial recovered ({nseeds}) than available susceptible ({len(i_susceptible)})"
                )
                i_recovered = np.random.choice(i_susceptible, size=nseeds, replace=False)
                self.model.people.state[i_recovered] = State.RECOVERED.value

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.R[tick + 1], self.model.people, State.RECOVERED, "Recovered")

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.R[tick + 1], self.model.people, State.RECOVERED, "Recovered")

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        return

    def plot(self):
        _fig, ax1 = plt.subplots()
        for node in range(self.model.nodes.count):
            ax1.plot(self.model.nodes.R[:, node], label=f"Node {node}")
        ax1.set_xlabel("Tick")
        ax1.set_ylabel("Recovered (by Node)")
        ax1.set_title("Recovered over Time by Node")
        ax1.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(np.sum(self.model.nodes.R, axis=1), color="black", linestyle="--", label="Total Recovered")
        ax2.set_ylabel("Total Recovered")
        ax2.legend(loc="upper right")
        plt.show()

        return


class RecoveredRS:
    """
    Recovered Component for SIRS/SEIRS Models (Waning Immunity)

    This component manages agents in the recovered state in models where immunity is temporary.
    It supports per-agent recovery timers, enabling individuals to return to the susceptible
    state after a configurable waning period. This is essential for SEIRS/SIRS model dynamics.

    Agents transition from Recovered back to Susceptible after the waning immunity period (rtimer).
    Tracks number of agents losing immunity each tick in `model.nodes.newly_waned`.

    Responsibilities:
    - Initializes agents in the `RECOVERED` state using `model.scenario.R`
    - Assigns `rtimer` values to track the duration of immunity
    - Decrements `rtimer` each tick; transitions agents to `SUSCEPTIBLE` when `rtimer == 0`
    - Updates patch-level counts:
        • `R[t, i]`: number of recovered individuals in node `i` at time `t`
        • `waned[t, i]`: number of agents who re-entered susceptibility on time step `t`

    Required Inputs:
    - `model.scenario.R`: initial number of recovered individuals per node
    - `wandurdist`: a function sampling the waning immunity duration distribution
    - `wandurmin`: minimum duration of immunity (default = 1 time step)

    Outputs:
    - `model.people.rtimer`: per-agent countdown to immunity expiration
    - `model.nodes.R`: recovered count per patch per timestep
    - `model.nodes.waned`: number of immunity losses per patch per tick

    Step Behavior:
    - Agents with `state == RECOVERED` decrement `rtimer`
    - When `rtimer == 0`, they return to `SUSCEPTIBLE`
    - `R` and `S` counts are updated to reflect this transition
    - `waned[t]` logs the number of agents who lost immunity on time step `t`

    Validation:
    - Ensures population conservation and consistency between agent states and patch totals
    - Detects unexpected changes in `R` or invalid transitions

    Plotting:
    The `plot()` method provides two views:
    1. Per-node recovered trajectories
    2. Total recovered and waned agents over time

    Example:
        model.components = [
            SIR.Susceptible(model),
            SEIRS.Infectious(model, infdurdist, wandurdist),
            Exposed(model, ...),
            RecoveredRS(model, wandurdist),
        ]
    """

    def __init__(self, model, wandurdist, wandurmin=1, validating: bool = False):
        self.model = model
        self.validating = validating

        if not hasattr(self.model.people, "rtimer"):
            self.model.people.add_scalar_property("rtimer", dtype=np.uint16)
        self.model.nodes.add_vector_property("R", model.params.nticks + 1, dtype=np.int32)
        self.model.nodes.add_vector_property("newly_waned", model.params.nticks + 1, dtype=np.int32)

        self.model.nodes.R[0] = self.model.scenario.R

        self.wandurdist = wandurdist
        self.wandurmin = wandurmin

        # convenience
        nodeids = self.model.people.nodeid
        states = self.model.people.state

        for node in range(self.model.nodes.count):
            nseeds = self.model.scenario.R[node]
            if nseeds > 0:
                i_susceptible = np.nonzero((nodeids == node) & (states == State.SUSCEPTIBLE.value))[0]
                assert nseeds <= len(i_susceptible), (
                    f"Node {node} has more initial recovered ({nseeds}) than available susceptible ({len(i_susceptible)})"
                )
                i_recovered = np.random.choice(i_susceptible, size=nseeds, replace=False)
                self.model.people.state[i_recovered] = State.RECOVERED.value
                samples = dists.sample_floats(self.wandurdist, np.zeros(nseeds, dtype=np.float32))
                samples = np.round(samples)
                samples = np.maximum(samples, self.wandurmin).astype(self.model.people.rtimer.dtype)
                self.model.people.rtimer[i_recovered] = samples
                assert np.all(self.model.people.rtimer[i_recovered] > 0), (
                    f"Recovered individuals should have rtimer > 0 ({self.model.people.rtimer[i_recovered].min()=})"
                )

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_flow_vs_census(self.model.nodes.R[tick + 1], self.model.people, State.RECOVERED, "Recovered")
        _check_timer_active(self.model.people.state, State.RECOVERED.value, self.model.people.rtimer, "Recovered", "rtimer")
        _check_state_timer_consistency(self.model.people.state, State.RECOVERED.value, self.model.people.rtimer, "Recovered", "rtimer")

        alive = self.model.people.state != State.DECEASED.value
        self.rtimers_one = (alive) & (self.model.people.rtimer == 1)

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_flow_vs_census(self.model.nodes.R[tick + 1], self.model.people, State.RECOVERED, "Recovered")
        _check_timer_active(self.model.people.state, State.RECOVERED.value, self.model.people.rtimer, "Recovered", "rtimer")
        _check_state_timer_consistency(self.model.people.state, State.RECOVERED.value, self.model.people.rtimer, "Recovered", "rtimer")

        assert np.all(self.model.people.state[self.rtimers_one] == State.SUSCEPTIBLE.value), (
            "Individuals with rtimer == 1 before should now be susceptible."
        )

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        newly_waned_by_node = np.zeros((nb.get_num_threads(), self.model.nodes.count), dtype=np.int32)
        nb_timer_update(
            self.model.people.state,
            State.RECOVERED.value,
            self.model.people.rtimer,
            State.SUSCEPTIBLE.value,
            newly_waned_by_node,
            self.model.people.nodeid,
        )
        newly_waned_by_node = newly_waned_by_node.sum(axis=0).astype(self.model.nodes.S.dtype)  # Sum over threads

        # state(t+1) = state(t) + ∆state(t)
        self.model.nodes.R[tick + 1] -= newly_waned_by_node
        self.model.nodes.S[tick + 1] += newly_waned_by_node
        # Record today's ∆
        self.model.nodes.newly_waned[tick] = newly_waned_by_node

        return

    def plot(self):
        # First plot: Recovered over Time by Node
        _fig, ax1 = plt.subplots()
        for node in range(self.model.nodes.count):
            ax1.plot(self.model.nodes.R[:, node], label=f"Node {node}")
        ax1.set_xlabel("Tick")
        ax1.set_ylabel("Recovered (by Node)")
        ax1.set_title("Recovered over Time by Node")
        ax1.legend(loc="upper left")
        plt.show()

        # Second plot: Total Recovered and Total Waned over Time
        _fig, ax2 = plt.subplots()
        total_recovered = np.sum(self.model.nodes.R, axis=1)
        total_waned = np.sum(self.model.nodes.newly_waned, axis=1)
        ax2.plot(total_recovered, color="green", linestyle="--", label="Total Recovered")
        ax2.plot(total_waned, color="purple", linestyle="-.", label="Total Waned")
        ax2.set_xlabel("Tick")
        ax2.set_ylabel("Count")
        ax2.set_title("Total Recovered and Total Waned Over Time")
        ax2.legend(loc="upper right")
        plt.show()

        return


class TransmissionSIx:
    """
    Transmission Component for SI-Style Models (S → I Only, No Recovery)

    This component simulates the transmission process in simple epidemic models where
    agents move from the `SUSCEPTIBLE` to `INFECTIOUS` state and remain infectious
    indefinitely. It computes the force of infection (FOI) for each patch and applies
    it stochastically to susceptible agents.

    Agents transition from Susceptible to Infectious based on force of infection.
    Tracks number of new infections each tick in `model.nodes.newly_infected`.

    Responsibilities:
    - Computes per-node force of infection (`λ`) at each tick:
        λ = β * (I / N), with spatial coupling via a migration matrix
    - Applies probabilistic infection to susceptible agents using `nb_transmission_step`
    - Updates per-node `S` and `I` counts accordingly
    - Tracks new infections (incidence) and FOI values per node and tick

    Required Inputs:
    - `model.nodes.I[t]`: number of infectious agents per node at tick `t`
    - `model.nodes.S[t]`: number of susceptible agents per node at tick `t`
    - `model.params.beta`: transmission rate (global)
    - `model.network`: matrix of spatial coupling between nodes

    Outputs:
    - `model.nodes.forces[t, i]`: force of infection in node `i` at tick `t`
    - `model.nodes.incidence[t, i]`: number of new infections in node `i` at tick `t`

    Step Behavior:
    - Computes FOI (`λ`) for each node
    - Applies inter-node infection pressure via `model.network`
    - Converts FOI into a Bernoulli probability using: `p = 1 - exp(-λ)`
    - Infects susceptible agents probabilistically
    - Updates state and records incidence

    Validation:
    - Ensures consistency between state transitions and incidence records
    - Checks conservation of population in `S` and `I` states
    - Validates `incidence[t] == I[t+1] - I[t]`

    Plotting:
    The `plot()` method displays the force of infection over time per node.

    Example:
        model.components = [
            SIR.Susceptible(model),
            TransmissionSIx(model),
            InfectiousSI(model),
        ]
    """

    def __init__(self, model, seasonality: Union[ValuesMap, np.ndarray] = None, validating: bool = False):
        """
        Transmission Component for SI-Style Models (S → I Only, No Recovery)

        Args:
            model (Model): The epidemic model instance.
            seasonality (Union[ValuesMap, np.ndarray], optional): Seasonality modifier for transmission rate.
                Can be a ValuesMap or a precomputed array. Defaults to None.
            validating (bool): Enable component-level validation. Defaults to False.
        """
        self.model = model
        self.validating = validating
        self.model.nodes.add_vector_property("forces", model.params.nticks + 1, dtype=np.float32)
        self.model.nodes.add_vector_property("newly_infected", model.params.nticks + 1, dtype=np.int32)

        self.seasonality = seasonality if seasonality is not None else ValuesMap.from_scalar(1.0, model.params.nticks, model.nodes.count)

        return

    @staticmethod
    @nb.njit(
        # (nb.int8[:], nb.uint16[:], nb.float32[:], nb.uint32[:, :]),
        nogil=True,
        parallel=True,
        cache=True,
    )
    def nb_transmission_step(states, nodeids, ft, newly_infected_by_node):
        for i in nb.prange(len(states)):
            if states[i] == State.SUSCEPTIBLE.value:
                # Check for infection
                draw = np.random.rand()
                nid = nodeids[i]
                if draw < ft[nid]:
                    states[i] = State.INFECTIOUS.value
                    newly_infected_by_node[nb.get_thread_id(), nid] += 1

        return

    def prevalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")

        I = self.model.nodes.I  # noqa: E741
        assert np.all(self.model.nodes.newly_infected[tick] == (I[tick + 1] - I[tick])), (
            "Incidence does not match change in Infectious counts"
        )

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        ft = self.model.nodes.forces[tick]

        N = _get_total_population(self, tick)

        ft[:] = self.model.params.beta * self.seasonality[tick] * self.model.nodes.I[tick] / N
        transfer = ft[:, None] * self.model.network
        ft += transfer.sum(axis=0)
        ft -= transfer.sum(axis=1)
        ft = -np.expm1(-ft)  # Convert to probability of infection

        newly_infected_by_node = np.zeros((nb.get_num_threads(), self.model.nodes.count), dtype=np.int32)
        self.nb_transmission_step(
            self.model.people.state,
            self.model.people.nodeid,
            ft,
            newly_infected_by_node,
        )
        newly_infected_by_node = newly_infected_by_node.sum(axis=0)  # Sum over threads

        # state(t+1) = state(t) + ∆state(t)
        self.model.nodes.S[tick + 1] -= newly_infected_by_node
        self.model.nodes.I[tick + 1] += newly_infected_by_node
        # Record today's ∆
        self.model.nodes.newly_infected[tick] = newly_infected_by_node

        return

    def plot(self):
        for node in range(self.model.nodes.count):
            plt.plot(self.model.nodes.forces[:-1, node], label=f"Node {node}")
        plt.xlabel("Tick")
        plt.ylabel("Force of Infection")
        plt.title("Force of Infection over Time by Node")
        plt.legend()
        plt.show()

        return


class TransmissionSI:
    """
    Transmission Component for SIS/SIR/SIRS Models (S → I with Duration)

    This component simulates the transition from `SUSCEPTIBLE` to `INFECTIOUS` in
    models where infectious individuals have a finite infection duration (`itimer`).
    It supports full spatial coupling and allows infection durations to vary by node
    and tick.

    Agents transition from Susceptible to Infectious based on force of infection.
    Sets newly infectious agents' infection timers (itimer) based on `infdurdist` and `infdurmin`.
    Tracks number of new infections each tick in `model.nodes.newly_infected`.

    Responsibilities:
    - Computes force of infection (FOI) `λ = β * (I / N)` per patch each tick
    - Applies optional spatial coupling via `model.network` (infection pressure transfer)
    - Converts FOI into Bernoulli probabilities using `p = 1 - exp(-λ)`
    - Infects susceptible agents stochastically, assigning per-agent `itimer`
    - Updates patch-level susceptible (`S`) and infectious (`I`) counts
    - Records number of new infections per tick in `model.nodes.incidence`

    Required Inputs:
    - `model.params.beta`: transmission rate (global)
    - `model.network`: [n x n] matrix of transmission coupling
    - `infdurdist(tick, node)`: callable sampling the infectious duration distribution
    - `model.people.itimer`: preallocated per-agent infection timer

    Outputs:
    - `model.nodes.forces[t, i]`: computed FOI in node `i` at time `t`
    - `model.nodes.incidence[t, i]`: new infections in node `i` on time step `t`

    Step Behavior:
    - Computes FOI (`λ`) for each node
    - Applies inter-node infection pressure via `model.network`
    - Converts FOI into a Bernoulli probability using: `p = 1 - exp(-λ)`
    - Infects susceptible agents probabilistically
    - Updates state and records incidence

    Validation:
    - Ensures consistency between incidence and change in `I`
    - Checks for correct state and population accounting before and after tick

    Plotting:
    The `plot()` method visualizes per-node FOI (`λ`) over simulation time.

    Example:
        model.components = [
            SIR.Susceptible(model),
            TransmissionSI(model, infdurdist),
            InfectiousIR(model, infdurdist),
        ]
    """

    def __init__(
        self,
        model,
        infdurdist: Callable[[int, int], float],
        infdurmin: int = 1,
        seasonality: Union[ValuesMap, np.ndarray] = None,
        validating: bool = False,
    ):
        """
        Initializes the TransmissionSI component.

        Args:
            model (Model): The epidemiological model instance.
            infdurdist (Callable[[int, int], float]): A function that returns the infectious duration for a given tick and node.
            infdurmin (int): Minimum infectious duration.
            seasonality (Union[ValuesMap, np.ndarray], optional): Seasonality modifier for transmission rate. Defaults to None.
            validating (bool): Enable component-level validation. Defaults to False.
        """
        self.model = model
        self.validating = validating
        self.model.nodes.add_vector_property("forces", model.params.nticks + 1, dtype=np.float32)
        self.model.nodes.add_vector_property("newly_infected", model.params.nticks + 1, dtype=np.int32)

        self.infdurdist = infdurdist
        self.infdurmin = infdurmin

        # Default is no temporal or spatial variation in transmission
        self.seasonality = seasonality if seasonality is not None else ValuesMap.from_scalar(1.0, model.params.nticks, model.nodes.count)

        return

    def prevalidate_step(self, tick: int) -> None:
        # Check ahead because I->R and R->S transitions may have happened meaning S[tick] and I[tick] are "out of date"
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")

        self.prv_inext = self.model.nodes.I[tick + 1].copy()

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_flow_vs_census(self.model.nodes.I[tick + 1], self.model.people, State.INFECTIOUS, "Infectious")

        Inext = self.model.nodes.I[tick + 1]
        assert np.all(self.model.nodes.newly_infected[tick] == (Inext - self.prv_inext)), (
            "Incidence does not match change in Infectious counts"
        )

        return

    @staticmethod
    @nb.njit(nogil=True, parallel=True)
    def nb_transmission_step(states, nodeids, ft, newly_infected_by_node, itimers, infdurdist, infdurmin, tick):
        for i in nb.prange(len(states)):
            if states[i] == State.SUSCEPTIBLE.value:
                # Check for infection
                draw = np.random.rand()
                nid = nodeids[i]
                if draw < ft[nid]:
                    states[i] = State.INFECTIOUS.value
                    itimers[i] = np.maximum(np.round(infdurdist(tick, nid)), infdurmin)  # Set the infection timer
                    newly_infected_by_node[nb.get_thread_id(), nid] += 1

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        ft = self.model.nodes.forces[tick]

        N = _get_total_population(self, tick)

        ft[:] = self.model.params.beta * self.seasonality[tick] * self.model.nodes.I[tick] / N
        transfer = ft[:, None] * self.model.network
        ft += transfer.sum(axis=0)
        ft -= transfer.sum(axis=1)
        ft = -np.expm1(-ft)  # Convert to probability of infection

        newly_infected_by_node = np.zeros((nb.get_num_threads(), self.model.nodes.count), dtype=np.int32)
        self.nb_transmission_step(
            self.model.people.state,
            self.model.people.nodeid,
            ft,
            newly_infected_by_node,
            self.model.people.itimer,
            self.infdurdist,
            self.infdurmin,
            tick,
        )
        newly_infected_by_node = newly_infected_by_node.sum(axis=0).astype(self.model.nodes.S.dtype)  # Sum over threads

        # state(t+1) = state(t) + ∆state(t)
        self.model.nodes.S[tick + 1] -= newly_infected_by_node
        self.model.nodes.I[tick + 1] += newly_infected_by_node
        # Record today's ∆
        self.model.nodes.newly_infected[tick] = newly_infected_by_node

        return

    def plot(self):
        for node in range(self.model.nodes.count):
            plt.plot(self.model.nodes.forces[:, node], label=f"Node {node}")
        plt.xlabel("Tick")
        plt.ylabel("Force of Infection")
        plt.title("Force of Infection over Time by Node")
        plt.legend()
        plt.show()

        return


class TransmissionSE:
    """
    Transmission component for an SEIR/SEIRS model with S -> E transition and incubation duration.

        This component simulates the transition from `SUSCEPTIBLE` to `EXPOSED` in models
        where infection includes an incubation period before agents become infectious.
        It handles stochastic exposure based on per-node force of infection (FOI), and
        assigns individual incubation timers to newly exposed agents.

        Agents transition from Susceptible to Exposed based on force of infection.
        Sets newly exposed agents' infection timers (etimer) based on `expdurdist` and `expdurmin`.
        Tracks number of new infections each tick in `model.nodes.newly_infected`.

        Responsibilities:
        - Computes force of infection `λ = β * (I / N)` at each tick per node
        - Adjusts FOI using `model.network` for inter-node transmission coupling. Required but can be nullified by filling with all zeros.
        - Applies FOI to susceptible agents to determine exposure
        - Assigns incubation durations (`etimer`) to each newly exposed agent
        - Updates node-level counts for `S` and `E` and logs daily incidence

        Required Inputs:
        - `model.params.beta`: global transmission rate
        - `model.network`: [n x n] matrix for FOI migration
        - `expdurdist(tick, node)`: callable that samples the exposure/incubation duration distribution
        - `expdurmin`: minimum incubation period (default = 1)

        Outputs:
        - `model.nodes.forces[t, i]`: computed FOI in node `i` at tick `t`
        - `model.nodes.incidence[t, i]`: new exposures per node per day
        - `model.people.etimer`: per-agent incubation countdown

        Step Behavior:
        - Computes FOI (`λ`) for each node
        - Optionally applies inter-node infection pressure via `model.network`
        - Converts FOI into a Bernoulli probability using: `p = 1 - exp(-λ)`
        - Infects susceptible agents probabilistically
        - Updates state and records incidence

        Validation:
        - Validates consistency between agent states and patch-level counts before and after tick
        - Confirms that `incidence[t] == E[t+1] - E[t]`

        Plotting:
        The `plot()` method shows per-node FOI (`λ`) trajectories over time.

        Example:
            model.components = [
                SIR.Susceptible(model),
                TransmissionSE(model, expdurdist),
                Exposed(model, ...),
                InfectiousIR(model, ...),
                Recovered(model),
            ]
    """

    def __init__(
        self,
        model,
        expdurdist: Callable[[int, int], float],
        expdurmin: int = 1,
        seasonality: Union[ValuesMap, np.ndarray] = None,
        validating: bool = False,
    ):
        """
        Initializes the TransmissionSE component.

        Args:
            model (Model): The epidemiological model instance.
            expdurdist (Callable[[int, int], float]): A function that returns the incubation duration for a given tick and node.
            expdurmin (int): Minimum incubation duration.
            seasonality (Union[ValuesMap, np.ndarray], optional): Seasonality modifier for transmission rate. Defaults to None.
            validating (bool): Enable component-level validation. Defaults to False.
        """
        self.model = model
        self.validating = validating
        self.model.nodes.add_vector_property("forces", model.params.nticks + 1, dtype=np.float32)
        self.model.nodes.add_vector_property("newly_infected", model.params.nticks + 1, dtype=np.int32)

        self.expdurdist = expdurdist
        self.expdurmin = expdurmin

        # Default is no temporal or spatial variation in transmission
        self.seasonality = seasonality if seasonality is not None else ValuesMap.from_scalar(1.0, model.params.nticks, model.nodes.count)

        return

    def prevalidate_step(self, tick: int) -> None:
        # Check ahead because E->I and R->S transitions may have happened meaning S[tick] and E[tick] are "out of date"
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_flow_vs_census(self.model.nodes.E[tick + 1], self.model.people, State.EXPOSED, "Exposed")

        self.prv_enext = self.model.nodes.E[tick + 1].copy()

        return

    def postvalidate_step(self, tick: int) -> None:
        _check_flow_vs_census(self.model.nodes.S[tick + 1], self.model.people, State.SUSCEPTIBLE, "Susceptible")
        _check_flow_vs_census(self.model.nodes.E[tick + 1], self.model.people, State.EXPOSED, "Exposed")

        Enext = self.model.nodes.E[tick + 1]
        assert np.all(self.model.nodes.newly_infected[tick] == (Enext - self.prv_enext)), (
            "Incidence does not match change in Exposed counts"
        )

        return

    @staticmethod
    @nb.njit(nogil=True, parallel=True)
    def nb_transmission_step(states, nodeids, ft, newly_infected_by_node, etimers, expdurdist, expdurmin, tick):
        for i in nb.prange(len(states)):
            if states[i] == State.SUSCEPTIBLE.value:
                # Check for infection
                draw = np.random.rand()
                nid = nodeids[i]
                if draw < ft[nid]:
                    states[i] = State.EXPOSED.value
                    etimers[i] = np.maximum(np.round(expdurdist(tick, nid)), expdurmin)  # Set the exposure timer
                    newly_infected_by_node[nb.get_thread_id(), nid] += 1

        return

    @validate(pre=prevalidate_step, post=postvalidate_step)
    def step(self, tick: int) -> None:
        ft = self.model.nodes.forces[tick]

        N = _get_total_population(self, tick)

        ft[:] = self.model.params.beta * self.seasonality[tick] * self.model.nodes.I[tick] / N
        transfer = ft[:, None] * self.model.network
        ft += transfer.sum(axis=0)
        ft -= transfer.sum(axis=1)
        ft = -np.expm1(-ft)  # Convert to probability of infection

        newly_infected_by_node = np.zeros((nb.get_num_threads(), self.model.nodes.count), dtype=np.int32)
        self.nb_transmission_step(
            self.model.people.state,
            self.model.people.nodeid,
            ft,
            newly_infected_by_node,
            self.model.people.etimer,
            self.expdurdist,
            self.expdurmin,
            tick,
        )
        newly_infected_by_node = newly_infected_by_node.sum(axis=0).astype(self.model.nodes.S.dtype)  # Sum over threads

        # state(t+1) = state(t) + ∆state(t)
        self.model.nodes.S[tick + 1] -= newly_infected_by_node
        self.model.nodes.E[tick + 1] += newly_infected_by_node
        # Record today's ∆
        self.model.nodes.newly_infected[tick] = newly_infected_by_node

        return

    def plot(self):
        for node in range(self.model.nodes.count):
            plt.plot(self.model.nodes.forces[:, node], label=f"Node {node}")
        plt.xlabel("Tick")
        plt.ylabel("Force of Infection")
        plt.title("Force of Infection over Time by Node")
        plt.legend()
        plt.show()

        return


def _get_total_population(component, tick):
    """
    Helper to compute total population N at a given tick across all states.

    States, by default in laser.generic.Model are ['S', 'E', 'I', 'R'], but can be customized per component.
    E.g., a component may only care about ['S', 'I'] states or a component may add additional states like 'V' for vaccinated.

    Checks to see if the component has a 'states' attribute; if not, it uses the model's states.
    Requires the component to have a model property with nodes LaserFrame that contain current state counts.

    Args:
        component: The model component containing the model and states.
        tick: The current time tick.

    Returns:
        np.ndarray: Total population, shape = (#nodes,) array at the given tick.
    """

    states = component.states if hasattr(component, "states") else component.model.states

    # Use susceptible property as basis for N shape and dtype
    # Assumes that model at least has a susceptible property "S"
    N = np.zeros_like(component.model.nodes.S[tick])

    for state in states:
        # Test for state, e.g., might not have 'E' or 'R' states
        if (prop := getattr(component.model.nodes, state, None)) is not None:
            N += prop[tick]

    return N


## Validation helper functions


def _check_flow_vs_census(flow, people, state, name):
    """Compare a given flow vector against the census counts by state."""
    assert np.all(flow == (_actual := np.bincount(people.nodeid, people.state == state.value, len(flow)))), (
        f"{name} census does not match {name} counts (by state)."
    )
    return


def _check_unchanged(previous, current, name):
    """Check that a given array is unchanged after a step."""
    assert np.all(current == previous), f"{name} should be unchanged after step()."
    return


def _check_timer_active(states, value, timers, state_name, timer_name):
    """Check that individuals in a given state have active (greater than zero) timers."""
    assert np.all(_test := (timers[states == value] > 0)), (
        f"{state_name} individuals should have {timer_name} > 0 ({timers[states == value].min()=})"
    )
    return


def _check_state_timer_consistency(states, value, timers, state_name, timer_name):
    """Check that only live individuals in a given state have active (greater than zero) timers."""
    alive = states != State.DECEASED.value
    active = timers > 0
    assert np.all(_test := (states[alive & active] == value)), f"Only {state_name} individuals should have {timer_name} > 0."
    return
