"""Benchmark Numba vs C++ vs Rust across three kernel scenarios and an
n_agents sweep.

For each scenario we run up to ten configurations:

    numba/parallel             numba/serial
    numba/parallel (xorshift)  numba/serial (xorshift)    (scenarios 2 + 3)
    c++/openmp                 c++/serial
    c++/openmp (xorshift)      c++/serial (xorshift)      (scenarios 2 + 3)
    rust/rayon                 rust/serial

The timed region includes:
  * the per-thread-bin counter allocation,
  * the kernel call,
  * for the parallel paths, the across-threads accumulation (counter.sum(axis=0)).

The per-iteration inputs (states/timers/nodeids/ft/itimers) are rebuilt
outside the timed region.

Run from the laser.generic root:
    bash tmp/perf/build.sh
    python tmp/perf/bench.py
"""

import ctypes
import os
import sys
import time
from pathlib import Path
from statistics import median

import numba as nb
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import numba_kernels as nk  # noqa: E402


# ---------------------------------------------------------------------------
# ctypes setup
# ---------------------------------------------------------------------------


i8_p = ctypes.POINTER(ctypes.c_int8)
u16_p = ctypes.POINTER(ctypes.c_uint16)
i32_p = ctypes.POINTER(ctypes.c_int32)
u64_p = ctypes.POINTER(ctypes.c_uint64)
f32_p = ctypes.POINTER(ctypes.c_float)
SAMPLER_T = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_int64, ctypes.c_int32)


def _load_lib(path: Path) -> ctypes.CDLL:
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `bash tmp/perf/build.sh` first."
        )
    return ctypes.CDLL(str(path))


libcpp = _load_lib(HERE / "libcpp_kernels.so")
librust = _load_lib(HERE / "librust_kernels.so")


def _declare(lib: ctypes.CDLL, name: str, argtypes, restype=None) -> None:
    fn = getattr(lib, name)
    fn.argtypes = argtypes
    fn.restype = restype


# Scenario 1 (no RNG): parallel + serial only.
_timer_update_parallel_argtypes = [
    i8_p, ctypes.c_int8, u16_p, ctypes.c_int8, i32_p, u16_p,
    ctypes.c_int64, ctypes.c_int64,
]
_timer_update_serial_argtypes = [
    i8_p, ctypes.c_int8, u16_p, ctypes.c_int8, i32_p, u16_p, ctypes.c_int64,
]
_declare(libcpp, "cpp_timer_update", _timer_update_parallel_argtypes)
_declare(libcpp, "cpp_timer_update_serial", _timer_update_serial_argtypes)
_declare(librust, "rust_timer_update", _timer_update_parallel_argtypes)
_declare(librust, "rust_timer_update_serial", _timer_update_serial_argtypes)

# Scenario 2: default RNG (mt19937 / SmallRng) + xorshift variants (C++ only).
_transmission_parallel_argtypes = [
    i8_p, u16_p, f32_p, i32_p, ctypes.c_int8, ctypes.c_int8,
    ctypes.c_int64, ctypes.c_int64,
]
_transmission_serial_argtypes = [
    i8_p, u16_p, f32_p, i32_p, ctypes.c_int8, ctypes.c_int8, ctypes.c_int64,
]
# Numba xs keeps the rng_state arg; C++ xs uses thread_local internally
# so the signature collapses to the default-RNG shape.
_transmission_xs_numba_argtypes_parallel = [
    i8_p, u16_p, f32_p, i32_p, u64_p, ctypes.c_int8, ctypes.c_int8,
    ctypes.c_int64, ctypes.c_int64,
]
_transmission_xs_numba_argtypes_serial = [
    i8_p, u16_p, f32_p, i32_p, u64_p, ctypes.c_int8, ctypes.c_int8,
    ctypes.c_int64,
]
# C++ xs argtypes are identical to the default-RNG ones.
_declare(libcpp, "cpp_transmission_step", _transmission_parallel_argtypes)
_declare(libcpp, "cpp_transmission_step_serial", _transmission_serial_argtypes)
_declare(libcpp, "cpp_transmission_step_xorshift", _transmission_parallel_argtypes)
_declare(
    libcpp, "cpp_transmission_step_xorshift_serial",
    _transmission_serial_argtypes,
)
_declare(librust, "rust_transmission_step", _transmission_parallel_argtypes)
_declare(librust, "rust_transmission_step_serial", _transmission_serial_argtypes)

# Scenario 3: same shape, with sampler.
_t3_parallel = [
    i8_p, u16_p, f32_p, i32_p, u16_p, SAMPLER_T,
    ctypes.c_int32, ctypes.c_int64, ctypes.c_int8, ctypes.c_int8,
    ctypes.c_int64, ctypes.c_int64,
]
_t3_serial = [
    i8_p, u16_p, f32_p, i32_p, u16_p, SAMPLER_T,
    ctypes.c_int32, ctypes.c_int64, ctypes.c_int8, ctypes.c_int8,
    ctypes.c_int64,
]
# Numba xs scenario-3 keeps rng_state; C++ xs is thread_local internally
# so its signature matches the default shape.
_t3_xs_numba_parallel = [
    i8_p, u16_p, f32_p, i32_p, u16_p, u64_p, SAMPLER_T,
    ctypes.c_int32, ctypes.c_int64, ctypes.c_int8, ctypes.c_int8,
    ctypes.c_int64, ctypes.c_int64,
]
_t3_xs_numba_serial = [
    i8_p, u16_p, f32_p, i32_p, u16_p, u64_p, SAMPLER_T,
    ctypes.c_int32, ctypes.c_int64, ctypes.c_int8, ctypes.c_int8,
    ctypes.c_int64,
]
_declare(libcpp, "cpp_transmission_step_with_sampler", _t3_parallel)
_declare(libcpp, "cpp_transmission_step_with_sampler_serial", _t3_serial)
_declare(libcpp, "cpp_transmission_step_with_sampler_xorshift", _t3_parallel)
_declare(libcpp, "cpp_transmission_step_with_sampler_xorshift_serial", _t3_serial)
_declare(librust, "rust_transmission_step_with_sampler", _t3_parallel)
_declare(librust, "rust_transmission_step_with_sampler_serial", _t3_serial)

# Sampler symbols + thread-count queries.
_declare(libcpp, "cpp_sampler", [ctypes.c_int64, ctypes.c_int32], ctypes.c_double)
_declare(librust, "rust_sampler", [ctypes.c_int64, ctypes.c_int32], ctypes.c_double)
_declare(libcpp, "cpp_num_threads", [], ctypes.c_int64)
_declare(librust, "rust_num_threads", [], ctypes.c_int64)

cpp_sampler_ptr = ctypes.cast(libcpp.cpp_sampler, SAMPLER_T)
rust_sampler_ptr = ctypes.cast(librust.rust_sampler, SAMPLER_T)

NUMBA_THREADS = nb.get_num_threads()
CPP_THREADS = int(libcpp.cpp_num_threads())
RUST_THREADS = int(librust.rust_num_threads())


# ---------------------------------------------------------------------------
# Data setup
# ---------------------------------------------------------------------------


def ptr(arr, ctype):
    return arr.ctypes.data_as(ctypes.POINTER(ctype))


def setup_data(n_agents: int, n_nodes: int, seed: int = 42) -> dict:
    """Build a fresh per-iteration input set. Same seed -> reproducible."""
    rng = np.random.default_rng(seed)
    return {
        "states": rng.integers(0, 2, size=n_agents, dtype=np.int8),
        "timers": rng.integers(1, 11, size=n_agents, dtype=np.uint16),
        "nodeids": rng.integers(0, n_nodes, size=n_agents, dtype=np.uint16),
        "ft": (rng.random(n_nodes) * 0.001).astype(np.float32),
        "itimers": np.zeros(n_agents, dtype=np.uint16),
    }


def make_rng_state(size: int, seed: int) -> np.ndarray:
    """Seed `size` independent uint64 streams. Non-zero so xorshift can step."""
    rng = np.random.default_rng(seed)
    state = rng.integers(1, (1 << 63) - 1, size=size, dtype=np.uint64)
    return state


# ---------------------------------------------------------------------------
# Bench loop
# ---------------------------------------------------------------------------


def bench(label: str, n_iter: int, prepare, run) -> tuple[float, float]:
    """Returns (median_ms, mean_transitions). Prints inline."""
    for _ in range(3):
        run(prepare())

    times = []
    counts = []
    for _ in range(n_iter):
        inputs = prepare()
        t0 = time.perf_counter()
        per_node = run(inputs)
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000.0)
        counts.append(int(per_node.sum()))

    med_ms = median(times)
    mean_n = sum(counts) / len(counts)
    print(
        f"    {label:>28s}  median {med_ms:8.3f} ms  "
        f"mean transitions {mean_n:>12,.0f}"
    )
    return med_ms, mean_n


# ---------------------------------------------------------------------------
# Scenarios — each returns {label: (median_ms, mean_transitions)}
# ---------------------------------------------------------------------------


def run_scenario_1(n_agents: int, n_nodes: int, n_iter: int) -> dict:
    test_state = np.int8(1)
    new_state = np.int8(2)

    def prepare():
        d = setup_data(n_agents, n_nodes)
        d["states"][:] = 0
        d["states"][: n_agents // 2] = test_state
        d["timers"][: n_agents // 2] = 1
        return d

    results = {}

    def run_numba_parallel(d):
        counter = np.zeros((NUMBA_THREADS, n_nodes), dtype=np.int32)
        nk.numba_timer_update(
            d["states"], test_state, d["timers"], new_state, counter, d["nodeids"]
        )
        return counter.sum(axis=0)

    results["numba/parallel"] = bench("numba/parallel", n_iter, prepare, run_numba_parallel)

    def run_numba_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        nk.numba_timer_update_serial(
            d["states"], test_state, d["timers"], new_state, counter, d["nodeids"]
        )
        return counter

    results["numba/serial"] = bench("numba/serial", n_iter, prepare, run_numba_serial)

    def run_cpp_parallel(d):
        counter = np.zeros((CPP_THREADS, n_nodes), dtype=np.int32)
        libcpp.cpp_timer_update(
            ptr(d["states"], ctypes.c_int8), test_state,
            ptr(d["timers"], ctypes.c_uint16), new_state,
            ptr(counter, ctypes.c_int32),
            ptr(d["nodeids"], ctypes.c_uint16),
            ctypes.c_int64(n_agents), ctypes.c_int64(n_nodes),
        )
        return counter.sum(axis=0)

    results["c++/openmp"] = bench("c++/openmp", n_iter, prepare, run_cpp_parallel)

    def run_cpp_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        libcpp.cpp_timer_update_serial(
            ptr(d["states"], ctypes.c_int8), test_state,
            ptr(d["timers"], ctypes.c_uint16), new_state,
            ptr(counter, ctypes.c_int32),
            ptr(d["nodeids"], ctypes.c_uint16),
            ctypes.c_int64(n_agents),
        )
        return counter

    results["c++/serial"] = bench("c++/serial", n_iter, prepare, run_cpp_serial)

    def run_rust_parallel(d):
        counter = np.zeros((RUST_THREADS, n_nodes), dtype=np.int32)
        librust.rust_timer_update(
            ptr(d["states"], ctypes.c_int8), test_state,
            ptr(d["timers"], ctypes.c_uint16), new_state,
            ptr(counter, ctypes.c_int32),
            ptr(d["nodeids"], ctypes.c_uint16),
            ctypes.c_int64(n_agents), ctypes.c_int64(n_nodes),
        )
        return counter.sum(axis=0)

    results["rust/rayon"] = bench("rust/rayon", n_iter, prepare, run_rust_parallel)

    def run_rust_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        librust.rust_timer_update_serial(
            ptr(d["states"], ctypes.c_int8), test_state,
            ptr(d["timers"], ctypes.c_uint16), new_state,
            ptr(counter, ctypes.c_int32),
            ptr(d["nodeids"], ctypes.c_uint16),
            ctypes.c_int64(n_agents),
        )
        return counter

    results["rust/serial"] = bench("rust/serial", n_iter, prepare, run_rust_serial)
    return results


def run_scenario_2(n_agents: int, n_nodes: int, n_iter: int) -> dict:
    susceptible = np.int8(0)
    infectious = np.int8(1)

    # Persistent xorshift RNG state per backend. Allocated once outside the
    # timed region; the kernel mutates it in place across calls (matching
    # the thread_local pattern used by the default-RNG variants).
    rng_state_numba_par = make_rng_state(NUMBA_THREADS * 16, seed=11)
    rng_state_numba_ser = make_rng_state(1, seed=12)

    def prepare():
        d = setup_data(n_agents, n_nodes)
        d["states"][:] = susceptible
        d["ft"] = np.full(n_nodes, 0.01, dtype=np.float32)
        return d

    results = {}

    def run_numba_parallel(d):
        counter = np.zeros((NUMBA_THREADS, n_nodes), dtype=np.int32)
        nk.numba_transmission_step(
            d["states"], d["nodeids"], d["ft"], counter, susceptible, infectious,
        )
        return counter.sum(axis=0)

    results["numba/parallel"] = bench("numba/parallel", n_iter, prepare, run_numba_parallel)

    def run_numba_parallel_xs(d):
        counter = np.zeros((NUMBA_THREADS, n_nodes), dtype=np.int32)
        nk.numba_transmission_step_xorshift(
            d["states"], d["nodeids"], d["ft"], counter, rng_state_numba_par,
            susceptible, infectious,
        )
        return counter.sum(axis=0)

    results["numba/parallel (xs)"] = bench(
        "numba/parallel (xorshift)", n_iter, prepare, run_numba_parallel_xs
    )

    def run_numba_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        nk.numba_transmission_step_serial(
            d["states"], d["nodeids"], d["ft"], counter, susceptible, infectious,
        )
        return counter

    results["numba/serial"] = bench("numba/serial", n_iter, prepare, run_numba_serial)

    def run_numba_serial_xs(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        nk.numba_transmission_step_xorshift_serial(
            d["states"], d["nodeids"], d["ft"], counter, rng_state_numba_ser,
            susceptible, infectious,
        )
        return counter

    results["numba/serial (xs)"] = bench(
        "numba/serial (xorshift)", n_iter, prepare, run_numba_serial_xs
    )

    def run_cpp_parallel(d):
        counter = np.zeros((CPP_THREADS, n_nodes), dtype=np.int32)
        libcpp.cpp_transmission_step(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            susceptible, infectious,
            ctypes.c_int64(n_agents), ctypes.c_int64(n_nodes),
        )
        return counter.sum(axis=0)

    results["c++/openmp"] = bench("c++/openmp", n_iter, prepare, run_cpp_parallel)

    def run_cpp_parallel_xs(d):
        counter = np.zeros((CPP_THREADS, n_nodes), dtype=np.int32)
        libcpp.cpp_transmission_step_xorshift(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            susceptible, infectious,
            ctypes.c_int64(n_agents), ctypes.c_int64(n_nodes),
        )
        return counter.sum(axis=0)

    results["c++/openmp (xs)"] = bench(
        "c++/openmp (xorshift)", n_iter, prepare, run_cpp_parallel_xs
    )

    def run_cpp_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        libcpp.cpp_transmission_step_serial(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            susceptible, infectious,
            ctypes.c_int64(n_agents),
        )
        return counter

    results["c++/serial"] = bench("c++/serial", n_iter, prepare, run_cpp_serial)

    def run_cpp_serial_xs(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        libcpp.cpp_transmission_step_xorshift_serial(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            susceptible, infectious,
            ctypes.c_int64(n_agents),
        )
        return counter

    results["c++/serial (xs)"] = bench(
        "c++/serial (xorshift)", n_iter, prepare, run_cpp_serial_xs
    )

    def run_rust_parallel(d):
        counter = np.zeros((RUST_THREADS, n_nodes), dtype=np.int32)
        librust.rust_transmission_step(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            susceptible, infectious,
            ctypes.c_int64(n_agents), ctypes.c_int64(n_nodes),
        )
        return counter.sum(axis=0)

    results["rust/rayon"] = bench("rust/rayon", n_iter, prepare, run_rust_parallel)

    def run_rust_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        librust.rust_transmission_step_serial(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            susceptible, infectious,
            ctypes.c_int64(n_agents),
        )
        return counter

    results["rust/serial"] = bench("rust/serial", n_iter, prepare, run_rust_serial)
    return results


def run_scenario_3(n_agents: int, n_nodes: int, n_iter: int) -> dict:
    susceptible = np.int8(0)
    infectious = np.int8(1)
    infdurmin = np.int32(1)
    tick = np.int64(0)

    rng_state_numba_par = make_rng_state(NUMBA_THREADS * 16, seed=21)
    rng_state_numba_ser = make_rng_state(1, seed=22)

    def prepare():
        d = setup_data(n_agents, n_nodes)
        d["states"][:] = susceptible
        d["ft"] = np.full(n_nodes, 0.01, dtype=np.float32)
        d["itimers"][:] = 0
        return d

    results = {}

    def run_numba_parallel(d):
        counter = np.zeros((NUMBA_THREADS, n_nodes), dtype=np.int32)
        nk.numba_transmission_step_with_sampler(
            d["states"], d["nodeids"], d["ft"], counter, d["itimers"],
            nk.numba_sampler, infdurmin, tick, susceptible, infectious,
        )
        return counter.sum(axis=0)

    results["numba/parallel"] = bench("numba/parallel", n_iter, prepare, run_numba_parallel)

    def run_numba_parallel_xs(d):
        counter = np.zeros((NUMBA_THREADS, n_nodes), dtype=np.int32)
        nk.numba_transmission_step_with_sampler_xorshift(
            d["states"], d["nodeids"], d["ft"], counter, d["itimers"],
            rng_state_numba_par, nk.numba_sampler, infdurmin, tick,
            susceptible, infectious,
        )
        return counter.sum(axis=0)

    results["numba/parallel (xs)"] = bench(
        "numba/parallel (xorshift)", n_iter, prepare, run_numba_parallel_xs
    )

    def run_numba_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        nk.numba_transmission_step_with_sampler_serial(
            d["states"], d["nodeids"], d["ft"], counter, d["itimers"],
            nk.numba_sampler, infdurmin, tick, susceptible, infectious,
        )
        return counter

    results["numba/serial"] = bench("numba/serial", n_iter, prepare, run_numba_serial)

    def run_numba_serial_xs(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        nk.numba_transmission_step_with_sampler_xorshift_serial(
            d["states"], d["nodeids"], d["ft"], counter, d["itimers"],
            rng_state_numba_ser, nk.numba_sampler, infdurmin, tick,
            susceptible, infectious,
        )
        return counter

    results["numba/serial (xs)"] = bench(
        "numba/serial (xorshift)", n_iter, prepare, run_numba_serial_xs
    )

    def run_cpp_parallel(d):
        counter = np.zeros((CPP_THREADS, n_nodes), dtype=np.int32)
        libcpp.cpp_transmission_step_with_sampler(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            ptr(d["itimers"], ctypes.c_uint16),
            cpp_sampler_ptr, infdurmin, tick,
            susceptible, infectious,
            ctypes.c_int64(n_agents), ctypes.c_int64(n_nodes),
        )
        return counter.sum(axis=0)

    results["c++/openmp"] = bench("c++/openmp", n_iter, prepare, run_cpp_parallel)

    def run_cpp_parallel_xs(d):
        counter = np.zeros((CPP_THREADS, n_nodes), dtype=np.int32)
        libcpp.cpp_transmission_step_with_sampler_xorshift(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            ptr(d["itimers"], ctypes.c_uint16),
            cpp_sampler_ptr, infdurmin, tick,
            susceptible, infectious,
            ctypes.c_int64(n_agents), ctypes.c_int64(n_nodes),
        )
        return counter.sum(axis=0)

    results["c++/openmp (xs)"] = bench(
        "c++/openmp (xorshift)", n_iter, prepare, run_cpp_parallel_xs
    )

    def run_cpp_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        libcpp.cpp_transmission_step_with_sampler_serial(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            ptr(d["itimers"], ctypes.c_uint16),
            cpp_sampler_ptr, infdurmin, tick,
            susceptible, infectious,
            ctypes.c_int64(n_agents),
        )
        return counter

    results["c++/serial"] = bench("c++/serial", n_iter, prepare, run_cpp_serial)

    def run_cpp_serial_xs(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        libcpp.cpp_transmission_step_with_sampler_xorshift_serial(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            ptr(d["itimers"], ctypes.c_uint16),
            cpp_sampler_ptr, infdurmin, tick,
            susceptible, infectious,
            ctypes.c_int64(n_agents),
        )
        return counter

    results["c++/serial (xs)"] = bench(
        "c++/serial (xorshift)", n_iter, prepare, run_cpp_serial_xs
    )

    def run_rust_parallel(d):
        counter = np.zeros((RUST_THREADS, n_nodes), dtype=np.int32)
        librust.rust_transmission_step_with_sampler(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            ptr(d["itimers"], ctypes.c_uint16),
            rust_sampler_ptr, infdurmin, tick,
            susceptible, infectious,
            ctypes.c_int64(n_agents), ctypes.c_int64(n_nodes),
        )
        return counter.sum(axis=0)

    results["rust/rayon"] = bench("rust/rayon", n_iter, prepare, run_rust_parallel)

    def run_rust_serial(d):
        counter = np.zeros(n_nodes, dtype=np.int32)
        librust.rust_transmission_step_with_sampler_serial(
            ptr(d["states"], ctypes.c_int8),
            ptr(d["nodeids"], ctypes.c_uint16),
            ptr(d["ft"], ctypes.c_float),
            ptr(counter, ctypes.c_int32),
            ptr(d["itimers"], ctypes.c_uint16),
            rust_sampler_ptr, infdurmin, tick,
            susceptible, infectious,
            ctypes.c_int64(n_agents),
        )
        return counter

    results["rust/serial"] = bench("rust/serial", n_iter, prepare, run_rust_serial)
    return results


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------


def _human_n(n: int) -> str:
    if n >= 1_000_000:
        return f"{n // 1_000_000:>3d}M"
    if n >= 1_000:
        return f"{n // 1_000:>3d}K"
    return f"{n:>4d}"


def print_summary(results: dict, sweep_ns: list) -> None:
    print("\n" + "=" * 96)
    print("SUMMARY (median ms)")
    print("=" * 96)

    for scenario_id in (1, 2, 3):
        if scenario_id == 1:
            title = "Scenario 1: timer_update"
        elif scenario_id == 2:
            title = "Scenario 2: transmission_step"
        else:
            title = "Scenario 3: transmission_step + sampler"
        print(f"\n{title}")

        # Stable label order from the first n_agents result.
        first_n = sweep_ns[0]
        labels = list(results[scenario_id][first_n].keys())

        header = f"  {'':>22s}"
        for n in sweep_ns:
            header += f"  {_human_n(n):>10s}"
        print(header)

        for label in labels:
            row = f"  {label:>22s}"
            for n in sweep_ns:
                med, _ = results[scenario_id][n][label]
                row += f"  {med:>8.2f} ms"
            print(row)


def main() -> None:
    sweep_ns = [100_000, 1_000_000, 10_000_000, 100_000_000]
    n_nodes = 100

    # Iteration counts roughly amortize to constant wall time per kernel.
    iter_for = {
        100_000: 200,
        1_000_000: 100,
        10_000_000: 25,
        100_000_000: 8,
    }

    print("Numba vs C++ vs Rust kernel comparison")
    print(
        f"Threads:  numba={NUMBA_THREADS}  c++/openmp={CPP_THREADS}  "
        f"rust/rayon={RUST_THREADS}"
    )
    print(
        f"Sweep:    n_agents = {', '.join(f'{n:,}' for n in sweep_ns)}, "
        f"n_nodes = {n_nodes}"
    )
    print("Timed:    counter alloc + kernel + (parallel) per-thread accumulation.")
    print()

    results: dict = {1: {}, 2: {}, 3: {}}

    for n_agents in sweep_ns:
        n_iter = iter_for[n_agents]
        for scenario_id, scenario_fn in (
            (1, run_scenario_1),
            (2, run_scenario_2),
            (3, run_scenario_3),
        ):
            print(
                f"\n--- scenario {scenario_id}, n_agents = {n_agents:,}, "
                f"n_iter = {n_iter} ---"
            )
            results[scenario_id][n_agents] = scenario_fn(n_agents, n_nodes, n_iter)

    print_summary(results, sweep_ns)


if __name__ == "__main__":
    main()
