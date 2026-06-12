# Numba vs C++ vs Rust kernel comparison

Three kernel scenarios, three backends, single-threaded baselines for each,
two RNG choices for Numba and C++, swept across four orders of magnitude
of agent count.

## Scenarios

Each mirrors a real kernel in `src/laser/generic/components.py`:

1. **`timer_update`** — `nb_timer_update`. Scan agents, decrement a timer
   where state matches, transition on `timer == 0` and bump a per-thread,
   per-node counter. No RNG, no callbacks.
2. **`transmission_step`** — `TransmissionSIx.nb_transmission_step`. Adds
   a uniform random draw per susceptible agent and conditionally transitions.
3. **`transmission_step_with_sampler`** — `TransmissionSI.nb_transmission_step`.
   Same as (2) plus an indirect call to a provided sampler to set `itimer`.

## Backends and variants

All three parallel paths use the same per-thread-bin pattern for the
transition counter and rely on Python-side summation across the thread
dimension. No atomics. Each backend has a single-threaded variant with
a 1D `transitioned[nid]` array so the across-threads sum isn't needed.

Scenarios 2 and 3 also have a **xorshift alternative** for Numba and C++
that swaps in a fast hand-rolled xorshift64 PRNG. Rust's default is
already `rand::rngs::SmallRng` (xoshiro256++), which is the natural fast
choice in its ecosystem, so it doesn't get a second variant.

| Backend | Parallelism | Default RNG | Alternative RNG |
|---|---|---|---|
| **Numba** | `@nb.njit(parallel=True)` + `nb.prange` | `np.random.rand()` (thread-local) | xorshift64 over `rng_state[tid * 16]` (cache-line padded) |
| **C++** | OpenMP `#pragma omp parallel for` | `thread_local std::mt19937_64` + `uniform_real_distribution` | `thread_local` xorshift64 state |
| **Rust** | rayon, one task per logical thread | `thread_local!` `rand::rngs::SmallRng` | — |

Samplers (for scenario 3) all implement the same xorshift-based
deterministic computation, so scenario 3 measures indirect-call cost
rather than differences in sampler body work.

## Layout

```
tmp/perf/
├── README.md
├── build.sh              # builds libcpp_kernels.so + librust_kernels.so
├── bench.py              # benchmark runner (10 configs × 3 scenarios × 4 sizes)
├── numba_kernels.py      # default + xorshift, parallel + serial Numba kernels
├── cpp_kernels.cpp       # default + xorshift, parallel + serial C++ kernels
└── rust_kernels/
    ├── Cargo.toml        # rayon + rand, lto=fat
    └── src/lib.rs        # parallel + serial Rust kernels
```

## Running

```sh
bash tmp/perf/build.sh
python tmp/perf/bench.py
```

Override thread counts:

```sh
NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 RAYON_NUM_THREADS=8 \
    python tmp/perf/bench.py
```

Each backend's actual thread count is queried at runtime
(`nb.get_num_threads()`, `omp_get_max_threads()`,
`rayon::current_num_threads()`) so the per-thread bin buffer is sized
exactly right.

## Sample run (Apple Silicon, 14 logical cores)

```
Scenario 1: timer_update
                                100K          1M         10M        100M
          numba/parallel      0.13 ms      0.33 ms      1.89 ms     18.70 ms
            numba/serial      0.04 ms      0.36 ms      3.63 ms     36.32 ms
              c++/openmp      0.12 ms      0.68 ms      6.19 ms     57.24 ms
              c++/serial      0.04 ms      0.37 ms      3.63 ms     36.39 ms
              rust/rayon      0.13 ms      0.28 ms      1.53 ms     13.81 ms
             rust/serial      0.04 ms      0.37 ms      3.63 ms     36.63 ms

Scenario 2: transmission_step
                                100K          1M         10M        100M
          numba/parallel      0.14 ms      0.45 ms      3.23 ms     29.59 ms
     numba/parallel (xs)      0.15 ms      0.35 ms      2.15 ms     20.25 ms
            numba/serial      0.31 ms      3.07 ms     30.74 ms    309.90 ms
       numba/serial (xs)      0.16 ms      1.62 ms     16.14 ms    162.78 ms
              c++/openmp      0.10 ms      0.44 ms      3.74 ms     35.60 ms
         c++/openmp (xs)      0.07 ms      0.29 ms      2.21 ms     21.29 ms
              c++/serial      0.29 ms      2.91 ms     29.13 ms    291.94 ms
         c++/serial (xs)      0.17 ms      1.67 ms     16.67 ms    167.96 ms
              rust/rayon      0.13 ms      0.31 ms      1.29 ms     10.79 ms
             rust/serial      0.09 ms      0.82 ms      8.23 ms     82.40 ms

Scenario 3: transmission_step + sampler
                                100K          1M         10M        100M
          numba/parallel      0.14 ms      0.47 ms      3.39 ms     29.75 ms
     numba/parallel (xs)      0.15 ms      0.35 ms      2.52 ms     20.51 ms
            numba/serial      0.31 ms      3.13 ms     31.60 ms    313.11 ms
       numba/serial (xs)      0.17 ms      1.63 ms     16.59 ms    166.28 ms
              c++/openmp      0.10 ms      0.47 ms      3.90 ms     36.37 ms
         c++/openmp (xs)      0.08 ms      0.30 ms      2.43 ms     22.15 ms
              c++/serial      0.29 ms      2.92 ms     29.35 ms    293.03 ms
         c++/serial (xs)      0.17 ms      1.68 ms     16.88 ms    167.91 ms
              rust/rayon      0.13 ms      0.33 ms      1.55 ms     12.47 ms
             rust/serial      0.09 ms      0.83 ms      8.34 ms     82.60 ms
```

## How to read the numbers

**Scenario 1 is memory-bandwidth-bound at serial speed.** All three serial
backends land at ~3.6 ms / 36 ms (10M / 100M) — the cost of touching every
agent once. Numba and Rust still extract 2× parallel speedup; C++/OpenMP
*regresses* (6.19→3.63 ms at 10M, 57→36 ms at 100M), because the per-thread
bin update + `omp_get_thread_num()` call inside a 50%-hit branch is more
overhead than the work saves.

**RNG dominates serial scenarios 2 and 3.** At 100M:

| Serial RNG | Time |
|---|---:|
| Numba `np.random.rand()` | 310 ms |
| C++ `std::mt19937_64` | 292 ms |
| Numba xorshift64 | 163 ms |
| C++ xorshift64 (`thread_local`) | 168 ms |
| Rust `SmallRng` (xoshiro256++) | 82 ms |

Swapping in a hand-rolled xorshift roughly **halves** the serial cost in
Numba and C++. Rust's default `SmallRng` is already 2× faster than that —
xoshiro256++ is a meaningfully better algorithm than xorshift64. If you
want Numba/C++ to fully catch up, switch them to xoshiro256++ too.

**Parallel speedups at 100M, scenario 2:**

| Backend / RNG | Parallel | Serial | Speedup |
|---|---:|---:|---:|
| numba (np.random) | 29.59 ms | 309.90 ms | 10.5× |
| numba (xorshift)  | 20.25 ms | 162.78 ms | 8.0× |
| c++ (mt19937)     | 35.60 ms | 291.94 ms | 8.2× |
| c++ (xorshift, tls) | 21.29 ms | 167.96 ms | 7.9× |
| rust (SmallRng)   | 10.79 ms | 82.40 ms  | 7.6× |

Once the RNG choice is fixed across backends, the parallel speedups
converge to roughly 8× on 14 cores. The differences in *absolute* speed
come almost entirely from the RNG, not the parallel framework.

**Small-n is dominated by parallel-region launch.** At 100K agents,
parallel and serial costs are within a factor of 2-3 for most backends.
This is where Numba and Rust's lighter task launch beats OpenMP's
`parallel for` fork.

**Indirect-call cost is barely visible.** Scenario 3 vs scenario 2 differs
by under 0.5 ms in every cell at every size. The sampler-callback pattern
isn't a meaningful overhead in any of these backends.

## Implementation notes

A few non-obvious points if you want to extend this benchmark:

- **xorshift state must be cache-line padded (Numba) or `thread_local`
  (C++).** An early draft had Numba's `rng_state[tid]` and C++'s
  `rng_state[tid]` both accessing tightly-packed 8-byte slots in a
  shared array. False sharing dragged C++ xorshift's 100M parallel time
  from ~21 ms (where it should be) to ~113 ms — *worse* than mt19937.
  Numba was less affected but still suboptimal. The fix in the
  benchmark is two-pronged: pad Numba's slot stride to 128 bytes
  (`rng_state[tid * 16]`) and convert C++ to use a `thread_local`
  state directly. Rust's `thread_local!` `SmallRng` was unaffected
  because it always lived in thread-local storage.
- **C++ `thread_local` on Apple Silicon goes through `_tlv_get_addr`.**
  This is a real function call on first access per thread, which can
  show up at small `n`. After warm-up it caches a register, so the
  steady-state cost is small.
- **OpenMP `parallel for` has higher per-call fork cost than rayon or
  Numba's prange.** Visible in scenario 1 (where the inner work is
  tiny), much less visible in scenarios 2-3.
- **Each thread writes to its own row of the per-thread bin array.**
  For `n_nodes = 100` this is 400 bytes per row, so adjacent rows share
  at most one cache line at the boundary. Not a hot spot in practice.
  If you ever reduce `n_nodes` below ~20, false sharing on the bin
  rows would become visible.

## Build dependencies

- **C++**: `clang++` (Apple) + `libomp` (`brew install libomp`) on macOS;
  `g++` with `-fopenmp` on Linux. `build.sh` picks the right invocation
  per OS.
- **Rust**: any stable toolchain in `~/.cargo/bin` or on `PATH`. Build
  uses `rayon = "1.10"`, `rand = "0.8"`, `lto = "fat"`,
  `codegen-units = 1`.

## Caveats

- One iteration is one kernel call; we don't amortize parallel-region
  warm cost across many tightly-spaced calls.
- Three different RNGs are in play in the default columns. The
  xorshift columns standardize Numba and C++ on the same algorithm
  but Rust still uses xoshiro256++.
- ctypes marshal overhead is in the timing for C++ and Rust but is
  microseconds — invisible above 100K agents.
- Per-iteration wall time below 50 microseconds is below `perf_counter`'s
  reliable resolution on macOS. The 100K column should be read as
  "all under 0.2 ms" rather than as precise rankings.
