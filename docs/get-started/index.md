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

## LASER SIR Example – Annotated

This example demonstrates how to construct and run a simple **SIR** model in LASER using the `laser-core` and `laser-generic` libraries. It features:

- One spatial node (`1x1` grid)
- Poisson-distributed infectious periods
- Correct S → I → R transitions
- CSV output and plotting

---

### 1. Import Dependencies

```python
import numpy as np
import pandas as pd

from laser.core import PropertySet
from laser.core.distributions import poisson

from laser.generic.model import Model
from laser.generic import SIR
```

---

### 2. Define Parameters

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
rng = np.random.default_rng(params.seed)
```

---

### 3. Define Scenario (Single Patch)

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

---

### 4. Build Model

Initialize the `Model` using the scenario and parameters. The `.people` frame is automatically constructed with internal state fields.

```python
model = Model(scenario, params)
people = model.people  # Auto-generated LaserFrame for agents
```

---

### 5. Infectious Duration Distribution

We define a **Numba-wrapped Poisson distribution** for the infectious period using LASER’s distribution API.

```python
infectious_duration = poisson(params.mean_infectious_period)
```

---

### 6. Attach Components

LASER models are built from **modular components**, each responsible for a specific part of the disease process. These components are called **once per timestep** in the order provided.

The standard **SIR progression** involves:

1. Tracking the number of **Susceptible** agents (S)
2. Modeling the **Transmission** (S → I)
3. Modeling **Infectiousness and Recovery** (I → R)
4. Tracking the **Recovered** population (R)

We attach the components in this exact order to ensure state updates and population counts are handled consistently.

---

#### `SIR.Susceptible(model)`

This component:
- Initializes agents' state to `SUSCEPTIBLE` (code 0)
- Records the number of susceptible agents per node at each timestep
- **Does not modify** state on its own—it simply keeps track

No parameters or distributions are required.

---

#### `SIR.TransmissionSI(model, infdurdist=...)`

This is the **S → I transition** component:
- Computes **force of infection**:
  $$
  \lambda = \beta \cdot \frac{I}{N}
  $$
- For each susceptible agent, performs a Bernoulli trial with probability \( p = 1 - e^{-\lambda} \)
- If infected, agent's state becomes `INFECTIOUS`, and they are assigned an **infection duration** drawn from `infdurdist`
- The timer is stored in `itimer` (infection timer)

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

This means newly infected agents will remain infectious for a random number of days drawn from a **Poisson distribution with mean 7**.

**Alternative distributions** available in `laser.core.distributions`:
- `exponential(scale)`
- `gamma(shape, scale)`
- `lognormal(mean, sigma)`
- `constant_int(value)`
- `custom` (with tick/node-dependent logic)

> You must use a **Numba-compatible function** with signature `(tick: int, node: int) → float/int`

---

#### `SIR.InfectiousIR(model, infdurdist=...)`

This is the **I → R transition** component:
- Decrements each agent's `itimer` (infection timer) every timestep
- When `itimer == 0`, agent moves to `RECOVERED` state (code 3)
- Records the number of infectious and recovered agents per node at each timestep

This component **must use the same distribution** (`infdurdist`) as `TransmissionSI`, because it expects that `itimer` was set there.

---

#### `SIR.Recovered(model)`

This component:
- Tracks the number of recovered agents
- Propagates `R[t+1] = R[t] + new_recoveries`
- Does **not** initiate any transitions or timers

No parameters are needed.

---

#### Full Component Setup

```python
model.components = [
    SIR.Susceptible(model),                         # Track S
    SIR.TransmissionSI(model, infdurdist=infectious_duration),  # S → I
    SIR.InfectiousIR(model, infdurdist=infectious_duration),    # I → R
    SIR.Recovered(model),                           # Track R
]
```

> **Order matters** — make sure Susceptible and Recovered components wrap the transition steps.

---

#### Optional Enhancements

- You can replace `SIR.InfectiousIR` with `SIR.InfectiousIRS` for **waning immunity** (SIRS model).
- You can use `SIR.TransmissionSE` and `SIR.Exposed` components for SEIR models.
- Add importation (`Infect_Random_Agents`) or demography (`Births`, `Deaths`) as additional components.

```python
from laser.generic.importation import Infect_Random_Agents
model.components.append(Infect_Random_Agents(model))
```

---
---

### 7. Run Simulation

Run the simulation for the configured number of timesteps.

```python
model.run()
```

---

### 8. Extract SIR Time Series

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

---

### 9. Save to CSV

Export the results to disk for downstream analysis or plotting.

```python
df.to_csv("sir_timeseries.csv", index=False)
print("Saved sir_timeseries.csv")
```

---

### 10. Plot Results

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

---

### Notes

- This example uses the `laser-generic` components for SIR modeling.
- Timers (itimer) are automatically assigned using the specified distribution.
- The `grid()` function defines the spatial layout but is used here with a single node.

#### Using AI

For internal IDM users, you can use a pre-built AI interface, [JENNER-GPT](https://chatgpt.com/g/g-67e6b80cd3e88191ae01e058f9df665e-jenner-ic) to create your simulations.

<!-- should add some example prompts -->


## Tutorials

The [Tutorials](../tutorials/index.md) section begins with example code to demonstrate setting up [simple SIR models](../tutorials/sir.md) and gradually adding complexity. For an interactive experience, begin with the first Jupyter notebook tutorial [SI model with no demographics](../tutorials/notebooks/01_SI_nobirths_logistic_growth.ipynb) for a very simple implementation of the LASER model.
