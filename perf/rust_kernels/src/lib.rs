//! Rust kernels for the perf comparison.
//!
//! Parallelized with Rayon. The pattern is one rayon task per logical
//! thread (`rayon::current_num_threads()`), each processing a contiguous
//! chunk of the agent range. The chunk index doubles as the "per-thread
//! bin" row for the `transitioned` / `newly_infected_by_node` arrays, so
//! we avoid atomics and Python sums across the chunk dimension after the
//! kernel returns.
//!
//! `SmallRng` is used per-chunk (thread-local on the rayon worker), seeded
//! deterministically from the chunk index so streams are independent.

use std::cell::RefCell;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

// -------------------------------------------------------------------------
// Raw-pointer wrapper to let rayon closures capture FFI pointers.
// We promise (via the FFI contract) that each task only touches indices
// in its own chunk, and that the `transitioned` / `newly_infected_by_node`
// rows are disjoint across tasks. With that invariant, sharing pointers
// across tasks is safe.
// -------------------------------------------------------------------------

#[derive(Copy, Clone)]
struct SendPtr<T>(*mut T);
unsafe impl<T> Send for SendPtr<T> {}
unsafe impl<T> Sync for SendPtr<T> {}

impl<T> SendPtr<T> {
    // Method on the wrapper (taking `self` by value) so closures never
    // borrow the inner raw pointer field directly. Without this, the
    // borrow checker tracks accesses to `.0` as `&*mut T`, which isn't
    // Sync regardless of the unsafe impls above.
    #[inline(always)]
    unsafe fn add(self, n: usize) -> *mut T {
        unsafe { self.0.add(n) }
    }
}

// Pre-thread RNG, lazily seeded from the chunk index passed in by the
// caller. Stored in a thread_local so successive kernel invocations on
// the same worker continue the same stream.
thread_local! {
    static RNG: RefCell<Option<SmallRng>> = RefCell::new(None);
}

fn with_rng<F: FnOnce(&mut SmallRng) -> R, R>(chunk_idx: usize, f: F) -> R {
    RNG.with(|cell| {
        let mut slot = cell.borrow_mut();
        if slot.is_none() {
            let seed = 0xdeadbeef_u64.wrapping_add(
                (chunk_idx as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)
            );
            *slot = Some(SmallRng::seed_from_u64(seed));
        }
        f(slot.as_mut().unwrap())
    })
}

// Compute the [start, end) bounds for the chunk this task owns.
fn chunk_bounds(chunk_idx: usize, n_chunks: usize, n: usize) -> (usize, usize) {
    let chunk_size = (n + n_chunks - 1) / n_chunks;
    let start = chunk_idx * chunk_size;
    let end = ((chunk_idx + 1) * chunk_size).min(n);
    (start, end)
}


// -------------------------------------------------------------------------
// Sampler for scenario 3 — same xorshift as the Numba and C++ versions.
// -------------------------------------------------------------------------

#[no_mangle]
pub extern "C" fn rust_sampler(tick: i64, nid: i32) -> f64 {
    let mut x: u64 = (tick as u64)
        .wrapping_mul(2_654_435_761)
        .wrapping_add((nid as u64).wrapping_mul(40_503));
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    (x % 1000) as f64 / 1000.0 * 5.0 + 2.0
}


// -------------------------------------------------------------------------
// Scenario 1: timer update
// -------------------------------------------------------------------------

/// # Safety
/// `states`, `timers`, `node_ids` must point to arrays of length `n`.
/// `transitioned` must point to an array of length `n_chunks * n_nodes`,
/// where `n_chunks == rayon::current_num_threads()`.
#[no_mangle]
pub unsafe extern "C" fn rust_timer_update(
    states: *mut i8,
    test_state: i8,
    timers: *mut u16,
    new_state: i8,
    transitioned: *mut i32,
    node_ids: *const u16,
    n: i64,
    n_nodes: i64,
) {
    let n = n as usize;
    let n_nodes = n_nodes as usize;
    let n_chunks = rayon::current_num_threads();

    let states = SendPtr(states);
    let timers = SendPtr(timers);
    let transitioned = SendPtr(transitioned);
    let node_ids = SendPtr(node_ids as *mut u16);

    (0..n_chunks).into_par_iter().for_each(|chunk_idx| {
        let (start, end) = chunk_bounds(chunk_idx, n_chunks, n);
        for i in start..end {
            unsafe {
                let s = *states.add(i);
                if s == test_state {
                    let t = (*timers.add(i)).wrapping_sub(1);
                    *timers.add(i) = t;
                    if t == 0 {
                        *states.add(i) = new_state;
                        let nid = *node_ids.add(i) as usize;
                        *transitioned.add(chunk_idx * n_nodes + nid) += 1;
                    }
                }
            }
        }
    });
}


// -------------------------------------------------------------------------
// Scenario 2: transmission (uniform draw, no callback)
// -------------------------------------------------------------------------

/// # Safety
/// Same pointer-length invariants as `rust_timer_update`.
#[no_mangle]
pub unsafe extern "C" fn rust_transmission_step(
    states: *mut i8,
    node_ids: *const u16,
    ft: *const f32,
    newly_infected_by_node: *mut i32,
    susceptible_state: i8,
    infectious_state: i8,
    n: i64,
    n_nodes: i64,
) {
    let n = n as usize;
    let n_nodes = n_nodes as usize;
    let n_chunks = rayon::current_num_threads();

    let states = SendPtr(states);
    let node_ids = SendPtr(node_ids as *mut u16);
    let ft = SendPtr(ft as *mut f32);
    let infected = SendPtr(newly_infected_by_node);

    (0..n_chunks).into_par_iter().for_each(|chunk_idx| {
        let (start, end) = chunk_bounds(chunk_idx, n_chunks, n);
        with_rng(chunk_idx, |rng| {
            for i in start..end {
                unsafe {
                    if *states.add(i) == susceptible_state {
                        let draw: f64 = rng.gen();
                        let nid = *node_ids.add(i) as usize;
                        let ft_val = *ft.add(nid) as f64;
                        if draw < ft_val {
                            *states.add(i) = infectious_state;
                            *infected.add(chunk_idx * n_nodes + nid) += 1;
                        }
                    }
                }
            }
        });
    });
}


// -------------------------------------------------------------------------
// Scenario 3: transmission + sampler callback
// -------------------------------------------------------------------------

pub type SamplerFn = extern "C" fn(i64, i32) -> f64;

/// # Safety
/// Same pointer-length invariants as the other kernels. `sampler` must be
/// a valid C function pointer with the declared signature.
#[no_mangle]
pub unsafe extern "C" fn rust_transmission_step_with_sampler(
    states: *mut i8,
    node_ids: *const u16,
    ft: *const f32,
    newly_infected_by_node: *mut i32,
    itimers: *mut u16,
    sampler: SamplerFn,
    infdurmin: i32,
    tick: i64,
    susceptible_state: i8,
    infectious_state: i8,
    n: i64,
    n_nodes: i64,
) {
    let n = n as usize;
    let n_nodes = n_nodes as usize;
    let n_chunks = rayon::current_num_threads();

    let states = SendPtr(states);
    let node_ids = SendPtr(node_ids as *mut u16);
    let ft = SendPtr(ft as *mut f32);
    let infected = SendPtr(newly_infected_by_node);
    let itimers = SendPtr(itimers);

    (0..n_chunks).into_par_iter().for_each(|chunk_idx| {
        let (start, end) = chunk_bounds(chunk_idx, n_chunks, n);
        with_rng(chunk_idx, |rng| {
            for i in start..end {
                unsafe {
                    if *states.add(i) == susceptible_state {
                        let draw: f64 = rng.gen();
                        let nid = *node_ids.add(i) as usize;
                        let ft_val = *ft.add(nid) as f64;
                        if draw < ft_val {
                            *states.add(i) = infectious_state;
                            let sampled = sampler(tick, nid as i32).round() as i32;
                            let timer = sampled.max(infdurmin) as u16;
                            *itimers.add(i) = timer;
                            *infected.add(chunk_idx * n_nodes + nid) += 1;
                        }
                    }
                }
            }
        });
    });
}


// =========================================================================
// Serial variants — plain loops, no rayon, 1D accumulator (n_nodes,).
// Raw pointers don't need the SendPtr wrapper here since there's no closure
// boundary. The same with_rng(0, ...) reuses the per-thread RNG, which is
// fine for the calling thread.
// =========================================================================

/// # Safety
/// Same pointer-length invariants as the parallel variant, but `transitioned`
/// must point to an array of length `n_nodes` (1D).
#[no_mangle]
pub unsafe extern "C" fn rust_timer_update_serial(
    states: *mut i8,
    test_state: i8,
    timers: *mut u16,
    new_state: i8,
    transitioned: *mut i32,
    node_ids: *const u16,
    n: i64,
) {
    let n = n as usize;
    for i in 0..n {
        unsafe {
            let s = *states.add(i);
            if s == test_state {
                let t = (*timers.add(i)).wrapping_sub(1);
                *timers.add(i) = t;
                if t == 0 {
                    *states.add(i) = new_state;
                    let nid = *node_ids.add(i) as usize;
                    *transitioned.add(nid) += 1;
                }
            }
        }
    }
}

/// # Safety
/// Same pointer-length invariants as the parallel variant; `newly_infected_by_node`
/// is 1D of length `n_nodes`.
#[no_mangle]
pub unsafe extern "C" fn rust_transmission_step_serial(
    states: *mut i8,
    node_ids: *const u16,
    ft: *const f32,
    newly_infected_by_node: *mut i32,
    susceptible_state: i8,
    infectious_state: i8,
    n: i64,
) {
    let n = n as usize;
    with_rng(0, |rng| {
        for i in 0..n {
            unsafe {
                if *states.add(i) == susceptible_state {
                    let draw: f64 = rng.gen();
                    let nid = *node_ids.add(i) as usize;
                    let ft_val = *ft.add(nid) as f64;
                    if draw < ft_val {
                        *states.add(i) = infectious_state;
                        *newly_infected_by_node.add(nid) += 1;
                    }
                }
            }
        }
    });
}

/// # Safety
/// Same as the parallel scenario-3 variant; 1D accumulator.
#[no_mangle]
pub unsafe extern "C" fn rust_transmission_step_with_sampler_serial(
    states: *mut i8,
    node_ids: *const u16,
    ft: *const f32,
    newly_infected_by_node: *mut i32,
    itimers: *mut u16,
    sampler: SamplerFn,
    infdurmin: i32,
    tick: i64,
    susceptible_state: i8,
    infectious_state: i8,
    n: i64,
) {
    let n = n as usize;
    with_rng(0, |rng| {
        for i in 0..n {
            unsafe {
                if *states.add(i) == susceptible_state {
                    let draw: f64 = rng.gen();
                    let nid = *node_ids.add(i) as usize;
                    let ft_val = *ft.add(nid) as f64;
                    if draw < ft_val {
                        *states.add(i) = infectious_state;
                        let sampled = sampler(tick, nid as i32).round() as i32;
                        let timer = sampled.max(infdurmin) as u16;
                        *itimers.add(i) = timer;
                        *newly_infected_by_node.add(nid) += 1;
                    }
                }
            }
        }
    });
}


// -------------------------------------------------------------------------
// Thread-count query so Python can size the per-thread bin array exactly
// to the number of chunks rayon will spawn.
// -------------------------------------------------------------------------

#[no_mangle]
pub extern "C" fn rust_num_threads() -> i64 {
    rayon::current_num_threads() as i64
}
