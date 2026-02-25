# Overview

[![Documentation Build Status](https://github.com/laser-base/laser-generic/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/laser-base/laser-generic/actions/workflows/pages/pages-build-deployment)
[![build](https://github.com/laser-base/laser-generic/actions/workflows/github-actions.yml/badge.svg?branch=main)](https://github.com/laser-base/laser-generic/actions/workflows/github-actions.yml)
[![Coverage Status](https://codecov.io/gh/laser-base/laser-generic/branch/main/graphs/badge.svg?branch=main)](https://app.codecov.io/github/laser-base/laser-generic)
[![PyPI Package latest release](https://img.shields.io/pypi/v/laser-generic.svg)](https://test.pypi.org/project/laser-generic)
[![PyPI Wheel](https://img.shields.io/pypi/wheel/laser-generic.svg)](https://test.pypi.org/project/laser-generic)
[![Supported versions](https://img.shields.io/pypi/pyversions/laser-generic.svg)](https://test.pypi.org/project/laser-generic)
[![Supported implementations](https://img.shields.io/pypi/implementation/laser-generic.svg)](https://test.pypi.org/project/laser-generic)
[![Commits since latest release](https://img.shields.io/github/commits-since/InstituteforDiseaseModeling/laser-generic/v1.0.0.svg)](https://github.com/InstituteforDiseaseModeling/laser-generic/compare/v1.0.0...main)

LASER (Lightweight Agent Spatial modeling for ERadication) is a framework for building agent-based infectious disease models with an emphasis on spatial modeling and efficient computation at scale.

[`laser-generic`](https://github.com/laser-base/laser-generic) builds on top of [`laser-core`](https://github.com/laser-base/laser-core), offering a set of ready-to-use, generic disease model components (e.g., SI, SIS, SIR dynamics, births, deaths, vaccination).

* Free software: MIT license

## Getting Started and Documentation

We recommend using the [LASER documentation](https://laser.idmod.org/laser-generic) to familiarize yourself with the LASER disease modeling framework. However, the instructions below may be sufficient for those who want to jump right in.

### Installation

We recommend using [`uv`](https://github.com/astral-sh/uv) for faster, more reliable installs:

```
uv pip install laser-generic
```

Alternatively, you can use regular `pip`:

```
pip install laser-generic
```

To install the latest in-development version:

```
pip install https://github.com/laser-core/laser-generic/archive/main.zip
```

#### Using `laser-generic`

`laser-generic` can be used in your code after importing it into your project:

```python
import laser.generic as lg

print(lg.__version__)
```

### Development

1. clone the `laser-generic` repository with

```bash
git clone https://github.com/laser-base/laser-generic.git
```

2. install [`uv`](https://github.com/astral-sh/uv?tab=readme-ov-file#installation) _in your system [Python], i.e., _before_ creating and activating a virtual environment

3. install `tox` as a tool in `uv` with the `tox-uv` plugin with

```bash
uv tool install tox --with tox-uv
```

4. change to the `laser-generic` directoryh with

```bash
cd laser-generic
```

5. create a virtual environment for development with

```bash
uv venv .venv
```

6. activate the virtual environment with

**Mac or Linux:**

```bash
source .venv/bin/activate
```

**Windows:**

```sh
.venv\bin\Activate
```

#### Building Code in Development

**Option 1: build "live" code** - Python scripts using `laser.generic` will import the code directly from the repository clone. Edits to the source code will take effect upon restarting the Python environment and importing `laser.generic`.

with `pip`:

```bash
pip install -e ".[dev]"
```

with `uv`:

```bash
uv pip install -e ".[dev]"
```

**Option 2: build an installable package** - this package must be installed into an environment to be used and changes to the source code on disk will not be picked up by consumers of `laser.generic` until the package is rebuilt and reinstalled. However this process mirrors using `laser.generic` as a dependency better than the "live code" option above.

**Option 2A: build with `pip`**

- install the `build` package: `python3 -m pip install build`
- build the `laser-generic` Python packge: `python3 -m build`
- find the wheel file, `.whl`, in the `dist` directory: `ls -l dist`

**Option 2B: build with `uv`**

- `uv build`
- find the wheel file, `.whl`, in the `dist` directory: `ls -l dist`

#### Running Tests

Now you can run tests in the `tests` directory or run the entire check+docs+test suite with ```tox```. Running ```tox``` will run several consistency checks, build documentation, run tests against the supported versions of Python, and create a code coverage report based on the test suite. Note that the first run of ```tox``` may take a few minutes (~5). Subsequent runs should be quicker depending on the speed of your machine and the test suite (~2 minutes). You can use ```tox``` to run tests against a single version of Python with, for example, ```tox -e py312```.

-----

## Disclaimer

The code in this repository was developed by IDM and other collaborators to support our joint research on flexible agent-based modeling. We've made it publicly available under the MIT License to provide others with a better understanding of our research and an opportunity to build upon it for their own work. We make no representations that the code works as intended or that we will provide support, address issues that are found, or accept pull requests. You are welcome to create your own fork and modify the code to suit your own modeling needs as permitted under the MIT License.
