# Get started

LASER is intended for modeling diseases with straightforward and predictable transmission dynamics, such as influenza or measles.

## Installation

See the Home page for [installation instructions](../index.md#installation).

## Usage

To use the project:

```python
import laser.generic
laser.generic.compute(...)
```

## Create and run a simulation

LASER is designed to be modular and flexible to accommodate a variety of modeling needs. The example below demonstrates how to construct and run a simple **SIR** model in LASER using the `laser-core` and `laser-generic` libraries.
It features:

- One spatial node (`1x1` grid)
- Poisson-distributed infectious periods
- Correct S → I → R transitions
- CSV output and plotting


### 1. Import dependencies

```python
import numpy as np
import pandas as pd

from laser.core import PropertySet
from laser.core.distributions import poisson

from laser.generic.model import Model
from laser.generic import SIR
```

### 2. Define the parameters

We configure simulation-wide parameters using a `PropertySet`, including:

- Simulation length (`nticks`)
- Infection rate (`beta`)
- Average infectious period
- Number of initial infections
- RNG seed

```python
params = PropertySet({
    "nticks": 160,
    "beta": 0.8,                    # Per-day infection rate
    "mean_infectious_period": 7.0,  # Average duration of infectiousness
    "initial_infected": 10,
    "seed": 123,
})
```

### 3. Define the scenario (single patch)

Always use the `grid()` utility to create a scenario so that it's compliant with expectations downstream in the Model class. Here we use it even to create a 1x1 spatial node ("patch") with 50,000 people. The population is then split into S, I, and R:

```python
from laser.core.utils import grid
scenario = grid(
    M=1,
    N=1,
    population_fn=lambda r, c: 50_000
)

scenario["I"] = params.initial_infected
scenario["S"] = scenario["population"] - params.initial_infected
scenario["R"] = 0
```

### 4. Build the model

Initialize the `Model` using the scenario and parameters. The `.people` frame is automatically constructed with internal state fields.

```python
model = Model(scenario, params)
people = model.people  # Auto-generated LaserFrame for agents
```

### 5. Configure the infectious duration distribution

We define a **Numba-wrapped Poisson distribution** for the infectious period using LASER’s distribution API.

```python
infectious_duration = poisson(params.mean_infectious_period)
```

### 6. Attach components

LASER models are built from **modular components**, each responsible for a specific part of the disease process. Components are executed **once per timestep**, in the **order they are attached** to the model.

A standard **SIR** model is constructed from four conceptual steps:

- Tracking the number of susceptible agents (S)
- Modeling transmission from susceptible to infectious agents (S -> I)
- Modeling infectiousness and recovery (I -> R)
- Tracking the recovered population (R)

Correct ordering matters: components that record state must wrap components that change state, otherwise population counts will be inconsistent.

In this example, we attach the components directly from `laser.generic.components`, rather than using the `SIR` convenience submodule.

```python
from laser.generic.components import (
    Susceptible,
    TransmissionSI,
    InfectiousIR,
    Recovered,
)

model.components = [
    Susceptible(model),
    TransmissionSI(model, infdurdist=infdist),
    InfectiousIR(model, infdurdist=infdist),
    Recovered(model),
]
```

#### `Susceptible(model)`

This component:

- Initializes agents' infection state to **SUSCEPTIBLE** (state code `0`)
- Records the number of susceptible agents per node at each timestep
- **Does not modify state transitions** on its own

No parameters or distributions are required.

This component exists purely to track and record the susceptible population.

#### `TransmissionSI(model, infdurdist=...)`

This component implements the **S -> I transition**.

For each timestep, it:

- Computes the **force of infection**:
    $$
    \lambda = \beta \cdot \frac{I}{N}
    $$
- For each susceptible agent, performs a Bernoulli trial with probability:
    $$
    p = 1 - e^{-\lambda}
    $$
- If infection occurs:

    - The agent’s state is set to **INFECTIOUS** - An infection duration is drawn from `infdurdist`
    - The duration is stored in the agent’s `itimer` property

The `infdurdist` argument must be a **Numba-compatible distribution function**, for example:

```python
from laser.core.distributions import poisson
infdist = poisson(mean_infectious_period)
```

This component is responsible only for new infections; recovery is handled separately.

#### `InfectiousIR(model, infdurdist=...)`

This component handles the **I -> R transition**.

It:

- Decrements each infectious agent’s `itimer` each timestep
- Transitions agents to **RECOVERED** when their timer reaches zero
- Updates node-level counts for infectious and recovered populations

This component **must use the same `infdurdist`** as `TransmissionSI`, because it relies on the infection timers set during transmission.

#### `Recovered(model)`

This component:

- Tracks the number of recovered agents per node
- Updates recovered counts over time
- Does **not** initiate any transitions or timers

No parameters are required.

#### Important note on ordering

The recommended order is:

1. `Susceptible`
2. `TransmissionSI`
3. `InfectiousIR`
4. `Recovered`

This ensures that:

- Population counts are recorded consistently
- State transitions occur before recovery is tallied
- The invariant `S + I + R = N` is preserved at each timestep

##### Parameterization

The keyword argument:

```python
infdurdist=infectious_duration
```

is a Numba-wrapped distribution function. In this example, we use:

```python
from laser.core.distributions import poisson
infectious_duration = poisson(7.0)
```

This means newly infected agents will remain infectious for a random number of days drawn from a Poisson distribution with mean 7.

Alternative distributions available in `laser.core.distributions`:

- `exponential(scale)`
- `gamma(shape, scale)`
- `lognormal(mean, sigma)`
- `constant_int(value)`
- `custom` (with tick/node-dependent logic)

!!! note

     You must use a **Numba-compatible function** with signature `(tick: int, node: int) → float/int`


!!! note

    Order matters: make sure Susceptible and Recovered components wrap the transition steps.


#### Optional Enhancements

- You can replace `InfectiousIR` with `InfectiousIRS` for **waning immunity** (SIRS model).
- You can use `TransmissionSE` and `Exposed` components for SEIR models.
- Add importation (`Infect_Random_Agents`) or demography (`Births`, `Deaths`) as additional components.

```python
from laser.generic.importation import Infect_Random_Agents
model.components.append(Infect_Random_Agents(model))
```

### 7. Run the simulation

Run the simulation for the configured number of timesteps.

```python
model.run()
```


### 8. Extract SIR time series

Extract patch-level S, I, R results as a Pandas DataFrame.

```python
df = pd.DataFrame({
    "time": np.arange(params.nticks + 1),
    "S": model.nodes.S[:, 0],
    "I": model.nodes.I[:, 0],
    "R": model.nodes.R[:, 0],
})

print(df.head())
print("Peak infectious:", df["I"].max())
```


### 9. Save to CSV

Export the results to disk for downstream analysis or plotting.

```python
df.to_csv("sir_timeseries.csv", index=False)
print("Saved sir_timeseries.csv")
```


### 10. Plot results

Plot the trajectory of S, I, and R over time using `matplotlib`.

```python
import matplotlib.pyplot as plt

plt.plot(df["time"], df["S"], label="S")
plt.plot(df["time"], df["I"], label="I")
plt.plot(df["time"], df["R"], label="R")

plt.xlabel("Time")
plt.ylabel("Population")
plt.legend()
plt.grid(True)
plt.title("LASER SIR Example (1 node)")
plt.show()
```


## Using AI

For internal IDM users, you can use a pre-built AI interface, nicknamed [JENNER-GPT](https://chatgpt.com/g/g-67e6b80cd3e88191ae01e058f9df665e-jenner-ic), to create your simulations or ask questions about LASER. It is designed to know everything about LASER and can not only answer your general questions about the system, but also provide working code for components or for entire runnable scripts.


<!-- should add some example prompts -->


## Tutorials

The [Tutorials](../tutorials/index.md) section begins with example code to demonstrate setting up [simple SIR models](../tutorials/sir.md) and gradually adding complexity. For an interactive experience, begin with the first Jupyter notebook tutorial [SI model with no demographics](../tutorials/notebooks/01_SI_nobirths_logistic_growth.ipynb) for a very simple implementation of the LASER model.
