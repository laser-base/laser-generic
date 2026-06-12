"""Numba-compiled kernels for the perf comparison.

All three kernels use ``@nb.njit(parallel=True, nogil=True)`` and the
per-thread-bin accumulation pattern from the real laser-generic components
(``transitioned[nb.get_thread_id(), node_ids[i]] += 1``). Python-side code
sums across the thread dimension after the kernel returns.
"""

import numba as nb
import numpy as np


# ---------------------------------------------------------------------------
# Scenario 1: nb_timer_update — pure integer work, no random, no callback
# ---------------------------------------------------------------------------


@nb.njit(parallel=True, nogil=True, cache=True)
def numba_timer_update(
    states, test_state, timers, new_state, transitioned, node_ids
):
    for i in nb.prange(len(states)):
        if states[i] == test_state:
            timers[i] -= 1
            if timers[i] == 0:
                states[i] = new_state
                transitioned[nb.get_thread_id(), node_ids[i]] += 1


# ---------------------------------------------------------------------------
# Scenario 2: nb_transmission_step (SIx-style) — adds a uniform random draw
# ---------------------------------------------------------------------------


@nb.njit(parallel=True, nogil=True, cache=True)
def numba_transmission_step(
    states, nodeids, ft, newly_infected_by_node, susceptible, infectious
):
    for i in nb.prange(len(states)):
        if states[i] == susceptible:
            draw = np.random.rand()
            nid = nodeids[i]
            if draw < ft[nid]:
                states[i] = infectious
                newly_infected_by_node[nb.get_thread_id(), nid] += 1


# ---------------------------------------------------------------------------
# Sampler for scenario 3
# ---------------------------------------------------------------------------
#
# Deterministic xorshift-based pseudo-sample. Same algorithm is implemented
# in the C++ and Rust versions so the indirect-call cost is what we measure,
# not differences in the sampler body. No ``inline="always"`` so Numba is
# more likely to treat the call as indirect — matching how the real
# ``TransmissionSI`` calls a user-provided ``infdurdist``.


@nb.njit(nogil=True, cache=True)
def numba_sampler(tick, nid):
    x = np.uint64(np.int64(tick) * 2654435761 + np.int64(nid) * 40503)
    x ^= x << np.uint64(13)
    x ^= x >> np.uint64(17)
    x ^= x << np.uint64(5)
    return float(x % 1000) / 1000.0 * 5.0 + 2.0


# ---------------------------------------------------------------------------
# Scenario 3: nb_transmission_step (SI-style) — random draw + sampler call
# ---------------------------------------------------------------------------


@nb.njit(parallel=True, nogil=True, cache=True)
def numba_transmission_step_with_sampler(
    states,
    nodeids,
    ft,
    newly_infected_by_node,
    itimers,
    sampler,
    infdurmin,
    tick,
    susceptible,
    infectious,
):
    for i in nb.prange(len(states)):
        if states[i] == susceptible:
            draw = np.random.rand()
            nid = nodeids[i]
            if draw < ft[nid]:
                states[i] = infectious
                itimers[i] = max(int(round(sampler(tick, nid))), infdurmin)
                newly_infected_by_node[nb.get_thread_id(), nid] += 1


# ---------------------------------------------------------------------------
# Serial variants — same loop bodies, no nb.prange, no thread dimension on
# the accumulator. Useful for measuring the parallelization speedup.
# ---------------------------------------------------------------------------


@nb.njit(nogil=True, cache=True)
def numba_timer_update_serial(
    states, test_state, timers, new_state, transitioned, node_ids
):
    for i in range(len(states)):
        if states[i] == test_state:
            timers[i] -= 1
            if timers[i] == 0:
                states[i] = new_state
                transitioned[node_ids[i]] += 1


@nb.njit(nogil=True, cache=True)
def numba_transmission_step_serial(
    states, nodeids, ft, newly_infected_by_node, susceptible, infectious
):
    for i in range(len(states)):
        if states[i] == susceptible:
            draw = np.random.rand()
            nid = nodeids[i]
            if draw < ft[nid]:
                states[i] = infectious
                newly_infected_by_node[nid] += 1


@nb.njit(nogil=True, cache=True)
def numba_transmission_step_with_sampler_serial(
    states,
    nodeids,
    ft,
    newly_infected_by_node,
    itimers,
    sampler,
    infdurmin,
    tick,
    susceptible,
    infectious,
):
    for i in range(len(states)):
        if states[i] == susceptible:
            draw = np.random.rand()
            nid = nodeids[i]
            if draw < ft[nid]:
                states[i] = infectious
                itimers[i] = max(int(round(sampler(tick, nid))), infdurmin)
                newly_infected_by_node[nid] += 1


# ---------------------------------------------------------------------------
# Xorshift-RNG alternatives (scenarios 2 and 3 only).
#
# Hand-rolled xorshift64 instead of np.random.rand(). Per-thread state is
# passed in by the caller as a uint64 array of length n_threads (parallel)
# or length 1 (serial). The kernel reads/writes its own slot so streams
# stay independent. The same algorithm is implemented in cpp_kernels.cpp.
#
#   x ^= x << 13;  x ^= x >> 7;  x ^= x << 17;
#   draw = (x >> 11) * (1 / 2^53)
# ---------------------------------------------------------------------------

_XS_13 = np.uint64(13)
_XS_7 = np.uint64(7)
_XS_17 = np.uint64(17)
_XS_11 = np.uint64(11)
_XS_INV = 1.0 / float(1 << 53)


@nb.njit(parallel=True, nogil=True, cache=True)
def numba_transmission_step_xorshift(
    states, nodeids, ft, newly_infected_by_node, rng_state,
    susceptible, infectious,
):
    for i in nb.prange(len(states)):
        if states[i] == susceptible:
            tid = nb.get_thread_id()
            x = rng_state[tid * 16]
            x ^= x << _XS_13
            x ^= x >> _XS_7
            x ^= x << _XS_17
            rng_state[tid * 16] = x
            draw = np.float64(x >> _XS_11) * _XS_INV
            nid = nodeids[i]
            if draw < ft[nid]:
                states[i] = infectious
                newly_infected_by_node[tid, nid] += 1


@nb.njit(nogil=True, cache=True)
def numba_transmission_step_xorshift_serial(
    states, nodeids, ft, newly_infected_by_node, rng_state,
    susceptible, infectious,
):
    x = rng_state[0]
    for i in range(len(states)):
        if states[i] == susceptible:
            x ^= x << _XS_13
            x ^= x >> _XS_7
            x ^= x << _XS_17
            draw = np.float64(x >> _XS_11) * _XS_INV
            nid = nodeids[i]
            if draw < ft[nid]:
                states[i] = infectious
                newly_infected_by_node[nid] += 1
    rng_state[0] = x


@nb.njit(parallel=True, nogil=True, cache=True)
def numba_transmission_step_with_sampler_xorshift(
    states, nodeids, ft, newly_infected_by_node, itimers, rng_state,
    sampler, infdurmin, tick, susceptible, infectious,
):
    for i in nb.prange(len(states)):
        if states[i] == susceptible:
            tid = nb.get_thread_id()
            x = rng_state[tid * 16]
            x ^= x << _XS_13
            x ^= x >> _XS_7
            x ^= x << _XS_17
            rng_state[tid * 16] = x
            draw = np.float64(x >> _XS_11) * _XS_INV
            nid = nodeids[i]
            if draw < ft[nid]:
                states[i] = infectious
                itimers[i] = max(int(round(sampler(tick, nid))), infdurmin)
                newly_infected_by_node[tid, nid] += 1


@nb.njit(nogil=True, cache=True)
def numba_transmission_step_with_sampler_xorshift_serial(
    states, nodeids, ft, newly_infected_by_node, itimers, rng_state,
    sampler, infdurmin, tick, susceptible, infectious,
):
    x = rng_state[0]
    for i in range(len(states)):
        if states[i] == susceptible:
            x ^= x << _XS_13
            x ^= x >> _XS_7
            x ^= x << _XS_17
            draw = np.float64(x >> _XS_11) * _XS_INV
            nid = nodeids[i]
            if draw < ft[nid]:
                states[i] = infectious
                itimers[i] = max(int(round(sampler(tick, nid))), infdurmin)
                newly_infected_by_node[nid] += 1
    rng_state[0] = x
