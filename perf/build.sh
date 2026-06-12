#!/usr/bin/env bash
# Build the C++ and Rust shared libraries for the perf benchmark.
# Run from the laser.generic root: bash tmp/perf/build.sh

set -euo pipefail

cd "$(dirname "$0")"

# ---------------------------------------------------------------------------
# C++ — clang++ + libomp (Homebrew on macOS) or g++ on Linux.
# ---------------------------------------------------------------------------

echo "==> Building libcpp_kernels.so"

case "$(uname -s)" in
    Darwin)
        LIBOMP_PREFIX="$(brew --prefix libomp)"
        clang++ -O3 -fPIC -shared -std=c++17 \
            -Xpreprocessor -fopenmp \
            -I"$LIBOMP_PREFIX/include" \
            -L"$LIBOMP_PREFIX/lib" -lomp \
            cpp_kernels.cpp -o libcpp_kernels.so
        ;;
    *)
        g++ -O3 -fPIC -shared -std=c++17 -fopenmp \
            cpp_kernels.cpp -o libcpp_kernels.so
        ;;
esac

# ---------------------------------------------------------------------------
# Rust — cargo build --release inside the rust_kernels crate.
# ---------------------------------------------------------------------------

echo "==> Building librust_kernels.so"

(cd rust_kernels && "$HOME/.cargo/bin/cargo" build --release)

# Cargo names it librust_kernels.dylib on macOS, librust_kernels.so on Linux.
# Copy whichever exists to a stable name next to the bench script.
if [ -f rust_kernels/target/release/librust_kernels.dylib ]; then
    cp rust_kernels/target/release/librust_kernels.dylib librust_kernels.so
elif [ -f rust_kernels/target/release/librust_kernels.so ]; then
    cp rust_kernels/target/release/librust_kernels.so librust_kernels.so
else
    echo "could not find compiled Rust library" >&2
    exit 1
fi

echo "==> Done."
ls -l libcpp_kernels.so librust_kernels.so
