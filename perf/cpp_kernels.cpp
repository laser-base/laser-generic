// C++ kernels for the perf comparison.
//
// Parallelized with OpenMP. Uses per-thread bins for the transition counter
// (transitioned[tid * n_nodes + nid]) so we don't need atomics — matches the
// Numba pattern. Python sums across the thread dimension after the kernel
// returns.
//
// Build (macOS, libomp from Homebrew):
//   clang++ -O3 -fPIC -shared -std=c++17 \
//       -Xpreprocessor -fopenmp -I"$(brew --prefix libomp)/include" \
//       -L"$(brew --prefix libomp)/lib" -lomp \
//       cpp_kernels.cpp -o libcpp_kernels.so
//
// Each thread has its own std::mt19937_64 RNG (thread_local) seeded from
// the OpenMP thread index so different workers produce independent streams.

#include <cstdint>
#include <cmath>
#include <random>
#include <omp.h>

namespace {

// Per-thread RNG. Lazily seeded on first access from each worker thread
// using its omp_get_thread_num() so streams are distinct.
inline std::mt19937_64& thread_rng() {
    static thread_local std::mt19937_64 rng(
        0xdeadbeefULL +
        static_cast<uint64_t>(omp_get_thread_num()) * 0x9e3779b97f4a7c15ULL
    );
    return rng;
}

inline double uniform_draw() {
    static thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(thread_rng());
}

}  // namespace


extern "C" {

// -------------------------------------------------------------------------
// Sampler for scenario 3 — deterministic xorshift mirroring the Numba and
// Rust versions. Exposed so the bench can pass a function pointer to it.
// -------------------------------------------------------------------------

double cpp_sampler(int64_t tick, int32_t nid) {
    uint64_t x = static_cast<uint64_t>(
        static_cast<int64_t>(tick) * 2654435761LL +
        static_cast<int64_t>(nid)  * 40503LL
    );
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return static_cast<double>(x % 1000) / 1000.0 * 5.0 + 2.0;
}

typedef double (*sampler_t)(int64_t, int32_t);


// -------------------------------------------------------------------------
// Scenario 1: timer update
// -------------------------------------------------------------------------

void cpp_timer_update(
    int8_t* states,
    int8_t test_state,
    uint16_t* timers,
    int8_t new_state,
    int32_t* transitioned,    // (n_threads, n_nodes) row-major
    const uint16_t* node_ids,
    int64_t n,
    int64_t n_nodes
) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == test_state) {
            timers[i] -= 1;
            if (timers[i] == 0) {
                states[i] = new_state;
                const int tid = omp_get_thread_num();
                transitioned[tid * n_nodes + node_ids[i]] += 1;
            }
        }
    }
}


// -------------------------------------------------------------------------
// Scenario 2: transmission (uniform draw, no callback)
// -------------------------------------------------------------------------

void cpp_transmission_step(
    int8_t* states,
    const uint16_t* node_ids,
    const float* ft,
    int32_t* newly_infected_by_node,   // (n_threads, n_nodes) row-major
    int8_t susceptible_state,
    int8_t infectious_state,
    int64_t n,
    int64_t n_nodes
) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == susceptible_state) {
            const double draw = uniform_draw();
            const uint16_t nid = node_ids[i];
            if (draw < static_cast<double>(ft[nid])) {
                states[i] = infectious_state;
                const int tid = omp_get_thread_num();
                newly_infected_by_node[tid * n_nodes + nid] += 1;
            }
        }
    }
}


// -------------------------------------------------------------------------
// Scenario 3: transmission + sampler callback
// -------------------------------------------------------------------------

void cpp_transmission_step_with_sampler(
    int8_t* states,
    const uint16_t* node_ids,
    const float* ft,
    int32_t* newly_infected_by_node,
    uint16_t* itimers,
    sampler_t sampler,
    int32_t infdurmin,
    int64_t tick,
    int8_t susceptible_state,
    int8_t infectious_state,
    int64_t n,
    int64_t n_nodes
) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == susceptible_state) {
            const double draw = uniform_draw();
            const uint16_t nid = node_ids[i];
            if (draw < static_cast<double>(ft[nid])) {
                states[i] = infectious_state;
                int32_t timer = static_cast<int32_t>(
                    std::lround(sampler(tick, static_cast<int32_t>(nid)))
                );
                if (timer < infdurmin) timer = infdurmin;
                itimers[i] = static_cast<uint16_t>(timer);
                const int tid = omp_get_thread_num();
                newly_infected_by_node[tid * n_nodes + nid] += 1;
            }
        }
    }
}


// -------------------------------------------------------------------------
// Serial variants — no OpenMP, 1D accumulator (n_nodes,).
// -------------------------------------------------------------------------

void cpp_timer_update_serial(
    int8_t* states,
    int8_t test_state,
    uint16_t* timers,
    int8_t new_state,
    int32_t* transitioned,        // (n_nodes,) 1D
    const uint16_t* node_ids,
    int64_t n
) {
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == test_state) {
            timers[i] -= 1;
            if (timers[i] == 0) {
                states[i] = new_state;
                transitioned[node_ids[i]] += 1;
            }
        }
    }
}

void cpp_transmission_step_serial(
    int8_t* states,
    const uint16_t* node_ids,
    const float* ft,
    int32_t* newly_infected_by_node,   // (n_nodes,) 1D
    int8_t susceptible_state,
    int8_t infectious_state,
    int64_t n
) {
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == susceptible_state) {
            const double draw = uniform_draw();
            const uint16_t nid = node_ids[i];
            if (draw < static_cast<double>(ft[nid])) {
                states[i] = infectious_state;
                newly_infected_by_node[nid] += 1;
            }
        }
    }
}

void cpp_transmission_step_with_sampler_serial(
    int8_t* states,
    const uint16_t* node_ids,
    const float* ft,
    int32_t* newly_infected_by_node,   // (n_nodes,) 1D
    uint16_t* itimers,
    sampler_t sampler,
    int32_t infdurmin,
    int64_t tick,
    int8_t susceptible_state,
    int8_t infectious_state,
    int64_t n
) {
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == susceptible_state) {
            const double draw = uniform_draw();
            const uint16_t nid = node_ids[i];
            if (draw < static_cast<double>(ft[nid])) {
                states[i] = infectious_state;
                int32_t timer = static_cast<int32_t>(
                    std::lround(sampler(tick, static_cast<int32_t>(nid)))
                );
                if (timer < infdurmin) timer = infdurmin;
                itimers[i] = static_cast<uint16_t>(timer);
                newly_infected_by_node[nid] += 1;
            }
        }
    }
}

// -------------------------------------------------------------------------
// Thread-count query so Python can size the per-thread bin array exactly
// to omp_get_max_threads() (== what the parallel-for will fork).
// -------------------------------------------------------------------------

int64_t cpp_num_threads() {
    return static_cast<int64_t>(omp_get_max_threads());
}

}  // extern "C"


// =========================================================================
// Xorshift-RNG alternatives (scenarios 2 and 3 only).
//
// Uses thread_local storage so each OS worker has its own state, matching
// how the default mt19937_64 kernel does it. This avoids the cache-line
// contention we'd get with a shared per-thread state array.
//
// Same xorshift64 algorithm as the Numba and Rust versions:
//   x ^= x << 13;  x ^= x >> 7;  x ^= x << 17;
//   draw = (x >> 11) * (1 / 2^53)
// =========================================================================

namespace {

inline double xorshift_uniform_tls() {
    static thread_local uint64_t state =
        0xdeadbeefULL +
        static_cast<uint64_t>(omp_get_thread_num()) * 0x9e3779b97f4a7c15ULL;
    uint64_t x = state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    state = x;
    return static_cast<double>(x >> 11) * (1.0 / static_cast<double>(1ULL << 53));
}

}  // namespace

extern "C" {

void cpp_transmission_step_xorshift(
    int8_t* states,
    const uint16_t* node_ids,
    const float* ft,
    int32_t* newly_infected_by_node,    // (n_threads, n_nodes) row-major
    int8_t susceptible_state,
    int8_t infectious_state,
    int64_t n,
    int64_t n_nodes
) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == susceptible_state) {
            const double draw = xorshift_uniform_tls();
            const uint16_t nid = node_ids[i];
            if (draw < static_cast<double>(ft[nid])) {
                states[i] = infectious_state;
                const int tid = omp_get_thread_num();
                newly_infected_by_node[tid * n_nodes + nid] += 1;
            }
        }
    }
}

void cpp_transmission_step_xorshift_serial(
    int8_t* states,
    const uint16_t* node_ids,
    const float* ft,
    int32_t* newly_infected_by_node,    // (n_nodes,)
    int8_t susceptible_state,
    int8_t infectious_state,
    int64_t n
) {
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == susceptible_state) {
            const double draw = xorshift_uniform_tls();
            const uint16_t nid = node_ids[i];
            if (draw < static_cast<double>(ft[nid])) {
                states[i] = infectious_state;
                newly_infected_by_node[nid] += 1;
            }
        }
    }
}

void cpp_transmission_step_with_sampler_xorshift(
    int8_t* states,
    const uint16_t* node_ids,
    const float* ft,
    int32_t* newly_infected_by_node,
    uint16_t* itimers,
    sampler_t sampler,
    int32_t infdurmin,
    int64_t tick,
    int8_t susceptible_state,
    int8_t infectious_state,
    int64_t n,
    int64_t n_nodes
) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == susceptible_state) {
            const double draw = xorshift_uniform_tls();
            const uint16_t nid = node_ids[i];
            if (draw < static_cast<double>(ft[nid])) {
                states[i] = infectious_state;
                int32_t timer = static_cast<int32_t>(
                    std::lround(sampler(tick, static_cast<int32_t>(nid)))
                );
                if (timer < infdurmin) timer = infdurmin;
                itimers[i] = static_cast<uint16_t>(timer);
                const int tid = omp_get_thread_num();
                newly_infected_by_node[tid * n_nodes + nid] += 1;
            }
        }
    }
}

void cpp_transmission_step_with_sampler_xorshift_serial(
    int8_t* states,
    const uint16_t* node_ids,
    const float* ft,
    int32_t* newly_infected_by_node,
    uint16_t* itimers,
    sampler_t sampler,
    int32_t infdurmin,
    int64_t tick,
    int8_t susceptible_state,
    int8_t infectious_state,
    int64_t n
) {
    for (int64_t i = 0; i < n; ++i) {
        if (states[i] == susceptible_state) {
            const double draw = xorshift_uniform_tls();
            const uint16_t nid = node_ids[i];
            if (draw < static_cast<double>(ft[nid])) {
                states[i] = infectious_state;
                int32_t timer = static_cast<int32_t>(
                    std::lround(sampler(tick, static_cast<int32_t>(nid)))
                );
                if (timer < infdurmin) timer = infdurmin;
                itimers[i] = static_cast<uint16_t>(timer);
                newly_infected_by_node[nid] += 1;
            }
        }
    }
}

}  // extern "C"
