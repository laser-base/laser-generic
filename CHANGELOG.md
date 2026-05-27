# Changelog

All notable changes to this project are documented here.

## [1.1.0] - 2026-03-24

Fixed divide-by-zero error on empty population nodes and improved internal code quality.

### Bug Fixes

- Fixed divide-by-zero error that occurred when processing nodes with zero population.

### Internal

- Added tests to demonstrate and validate correct behavior with zero-population nodes.
- Updated tests to correctly import from utils and age_at_infection modules.
- Localized and cleaned up usage of np.maximum() for improved code clarity.
- Resolved linter warnings across the codebase.

## [1.0.1] - 2026-02-26

Maintenance release with expanded Python version support, dependency relaxation, and documentation improvements.

### Features

- Added support for Python 3.14 following Numba compatibility.
- Added Python 3.13 builds to the CI/CD pipeline.

### Bug Fixes

- Fixed magic number usage; standardized terminology to 'people' instead of 'population' throughout the codebase.
- Fixed pyproject.toml configuration for compatibility with released laser-core 1.0.1.

### Improvements

- Relaxed laser-core dependency constraint to allow any compatible 1.x release (previously pinned to 1.0.x).
- Suppressed non-critical matplotlib/pyparsing deprecation warnings in Python 3.9 test runs.
- Removed macOS x86 builds from CI/CD due to llvmlite installation incompatibility on GitHub runners.

### Documentation

- Added automatic installation of docs/requirements during documentation builds.
- Disabled automatic execution of Jupyter notebooks during doc builds to improve build times.
- Fixed formatting for images and equations in notebooks.
- Updated README with corrected GitHub Actions badge and other improvements.

## [1.0.0] - 2025-12-17

First stable release of the package, promoting from 0.1.0 to 1.0.0 with streamlined build and release infrastructure.

### Internal

- Simplified CI/CD pipeline to single-platform build, reflecting pure Python package requirements.
- Updated bump-my-version configuration to manage version numbers independently of setup.py.
- Prepared package for publication to PyPI.

## [0.1.0] - 2025-12-17

Initial public release of laser-generic v0.1.0, introducing a full suite of modular epidemiological components, comprehensive documentation, and a robust test suite.

### Features

- Added modular epidemiological model components: SI, SIR, SIRS, SEIR, SEIRS, and SIR with carrier state (SIRs).
- Added TransmissionSIx component with spatial and temporal seasonality support.
- Added BirthsByCBR component for birth modeling using crude birth rates.
- Added MortalityByCDR component for mortality modeling using crude death rates.
- Added MortalityByEstimator component for estimator-based mortality modeling.
- Added ConstantPopVitalDynamics component for constant-population vital dynamics.
- Added RoutineImmunizationEx component for routine immunization with extended configuration.
- Added campaign immunization and routine immunization classes.
- Added ValuesMap utility class with indexing support, scalar initialization, and time-series replication across multiple years.
- Added get_default_params() utility function for retrieving default simulation parameters.
- Added component-level validation flags and validation logic to model components.
- Added support for specifying which node states contribute to total simulation population in Model.
- Added Model propagation of current state counts on each tick to reduce ordering dependencies between components.
- Added temporal and spatial seasonality support with accompanying tests and notebooks.
- Added England and Wales measles dataset to the repository.
- Added flexible infection importation functions with start-time parameter support.
- Made contextily (base map provider) an optional dependency.
- Added script to convert Jupyter notebooks to Python files.
- Added support for JENNER-GPT description and capabilities in documentation.

### Breaking Changes

- Refactored package namespace from laser_generic to laser.generic across all source code, tests, and notebooks.
- Moved State from laser.generic.shared to laser.generic top-level namespace.
- Renamed nsteps to nticks in ValuesMap for consistency.
- Swapped nnodes and nticks argument order in ValuesMap.from_scalar() to match array indexing convention ([tick, node]).
- Removed VitalDynamicsNNN components in favor of dedicated births and mortality components.
- Removed grid() and estimate_capacity() in favor of laser.core implementations.
- Removed cache=True from Numba-compiled functions that accept function arguments.
- Renamed TransmissionSIX to TransmissionSIx to correctly reflect that X is not a model state.
- Changed flow accounting to use newly_xxx naming convention for all flow properties.
- Promoted laser.generic.models to the laser.generic top-level namespace.
- Updated laser.core dependency to 1.0.0.
- Moved new utility functions from newutils.py into utils.py; removed newutils.py.

### Bug Fixes

- Fixed alignment issue between contagion vector and network matrix in transmission step.
- Fixed initialization bug for infectious timers when using random-access initializer.
- Fixed itimer type to uint16 by default to prevent integer overflow.
- Fixed census functions to be called at the beginning of each tick, resolving vector alignment issues with S, I, R outputs.
- Fixed uint32/64 casting issues in census components for Python 3.11 and 3.12 compatibility.
- Fixed ri_exploration notebook bug.
- Fixed chi-squared test inputs in statistical tests.
- Fixed exception handling when laser.core raises an error while adding newborns.
- Fixed link to 'Create and run a simulation' section in documentation.
- Fixed broken links in CONTRIBUTING.md and other documentation pages.
- Fixed bulleted list formatting in CONTRIBUTING.md.
- Fixed Google Analytics tracking key name in mkdocs.yml.
- Corrected SortedQueue typo in README.

### Documentation

- Migrated documentation from Sphinx to MkDocs with Material theme.
- Added combined documentation build covering both laser-core and laser-generic.
- Added automated API reference generation using gen-files and literate-nav MkDocs plugins.
- Added architecture document describing system design and component model.
- Added guidelines for contributing and a code of conduct.
- Added comprehensive docstrings to all epidemiological component classes covering responsibilities, inputs, outputs, validation, step behavior, and usage examples.
- Added tutorial notebooks for seasonality, vital dynamics, routine immunization, and mortality.

### Testing

- Added comprehensive unit tests for SI, SIR, SIRS, SEIR, and SEIRS model components.
- Added tests for BirthsByCBR, MortalityByCDR, MortalityByEstimator, and RoutineImmunizationEx components.
- Added tests for ValuesMap, sample_dobs(), sample_dods(), and other utility functions.
- Added spatial and temporal seasonality tests.
- Added component validation tests.
- Added Kermack-McKendrick fitting tests for SIR model.
- Added beta sweep tests fitting simulated case traces against input parameters.
- Added BirthByCBR test verifying dobs tracking when track=True.
- Switched mortality test statistical method from KS to chi-squared.
- Loosened mortality test bounds from 1% to 3% to reduce flakiness.
- Added Python 3.14 to tox.ini test matrix.
- Updated macOS runner in GitHub Actions to use an x86-equipped machine.
- Marked known-failing model tests with xfail decorators.

### Internal

- Reorganized package source layout from src/laser_generic/ to src/laser/generic/ for proper namespace packaging.
- Switched to explicit imports throughout the codebase.
- Consolidated Numba-compiled update functions and made parallel accumulation functions thread-safe.
- Moved notebook data files into a dedicated data/ directory.
- Updated laser-core dependency through versions up to 1.0.0.
- Added numpy as an explicit direct dependency.
- Added Python 3.13 support note; Numba not yet supported on 3.13 so documentation builds use Python 3.12.
- Configured ruff for linting; excluded in-development notebooks from lint checks.
- Added prefer 'prng_seed' key in Model while maintaining backward compatibility with 'prngseed' and 'seed'.
- Made code of conduct link in CONTRIBUTING.md absolute for compatibility in both repository root and documentation contexts.
- Added mermaid.js and additional MkDocs extensions to docs/requirements.txt.
- Ensured all text files end with a newline per POSIX standard.
