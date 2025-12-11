# Census-Based SIR Simulation in LASER

This script implements a clean, componentized Susceptible-Infectious-Recovered (SIR) model using **census-based reporting**. It avoids all flow-based tracking (`new_infections`, `incidence`, `force_of_infection`) and instead updates population states directly, recording per-timestep counts from agent states.

---

## Imports and Configuration

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from laser.core import PropertySet
from laser.core.utils import grid
from laser.core.distributions import poisson
from laser.generic.model import Model
from laser.generic import SIR
from laser.generic.shared import State
from laser.core import distributions as dists
```

- `laser.generic.shared.State` defines agent-level disease states.
- `poisson(...)` returns a compiled duration sampler.
- We reuse `SIR.Susceptible` only to initialize state and `nodeid`.

---

## Parameters

```python
params = PropertySet({
    "nticks": 160,
    "beta": 0.8,
    "mean_infectious_period": 7.0,
    "initial_infected": 10,
    "seed": 123
})
```

Simulation parameters are wrapped in a `PropertySet` for consistency with LASER APIs.

---

## Scenario and Model Initialization

```python
scenario = grid(M=1, N=1, population_fn=lambda r, c: 50000)
scenario["S"] = scenario["population"] - params.initial_infected
scenario["R"] = 0

model = Model(scenario, params)
people = model.people
```

- We use a 1x1 grid with 50,000 agents in a single node.
- Initially susceptible count is offset by `initial_infected`, which we will seed manually.

---

## Add Node-Level Census Properties

```python
nnodes = model.nodes.count
nticks = params.nticks + 1  # include t=0

model.nodes.add_vector_property("I", length=nticks, dtype=np.int32)
model.nodes.add_vector_property("R", length=nticks, dtype=np.int32)
```

- LASER does not add `nodes.S`, `nodes.I`, or `nodes.R` by default unless using flow-based components.
- These vectors will store counts by timestep and node.

---

## Infectious Duration Distribution

```python
infdurdist = poisson(params.mean_infectious_period)
```

Returns a Numba-compiled function that generates random infectious durations.

---

## Component: `TransmissionCensusSI`

```python
class TransmissionCensusSI:
    ...
```

- Computes force of infection as `lambda_t = beta * I / N`.
- Infects susceptible agents probabilistically.
- Assigns state `S → I` and initializes `infection_timer`.

---

## Component: `InfectiousCensusIR`

```python
class InfectiousCensusIR:
    ...
```

- Decrements `infection_timer` for infectious agents.
- Transitions `I → R` once timers expire.

---

## Component: `CensusTracker`

```python
class CensusTracker:
    ...
```

- Counts `S`, `I`, and `R` per node per timestep.
- Updates `model.nodes.{S,I,R}[tick, node]` based on agent `state` and `nodeid`.

---

## Component: `SeedInitialInfections`

```python
class SeedInitialInfections:
    ...
```

This component is **defined but not used** — infection seeding is done manually near the end of the script.

---

## Component Registration

```python
model.components = [
    SIR.Susceptible(model),
    TransmissionCensusSI(model, beta=params.beta),
    InfectiousCensusIR(model),
    CensusTracker(model)
]
```

- We reuse `SIR.Susceptible` to create `people.state` and `nodeid`.
- All other components are custom-built to avoid flow-based tracking.

---

## Manual Infection Seeding

```python
infected_ix = rng.choice(people.count, size=initial, replace=False)
people.state[infected_ix] = State.INFECTIOUS.value
```

- Directly assigns a few agents to `INFECTIOUS` before the simulation starts.
- Could be moved to `SeedInitialInfections`, but left manual for clarity.

---

## Run Simulation

```python
model.run()
```

LASER runs `component.step(tick)` for each tick `0..nticks-1`.

---

## Save Results to CSV

```python
df = pd.DataFrame({
    "time": np.arange(params.nticks),
    "S": model.nodes.S[:params.nticks, 0],
    "I": model.nodes.I[:params.nticks, 0],
    "R": model.nodes.R[:params.nticks, 0],
})
df.to_csv("sir_census.csv", index=False)
```

- Slices only ticks `0..159` (no extra row).
- Keeps output in sync with actual simulation ticks.

---

## Plot Results

```python
plt.plot(df["time"], df["S"], label="S")
plt.plot(df["time"], df["I"], label="I")
plt.plot(df["time"], df["R"], label="R")
...
plt.show()
```

- Classic SIR curve using matplotlib.
- Labels and grid added for readability.

---

## Design Notes

- **Census-based modeling** avoids all flows (`new_infections`, `incidence`, etc.).
- **Agent state is the only source of truth**; node-level stats are derived by counting.
- **Recovery and transmission separated** for clarity and testing.
- **SIR.Susceptible** reused only to avoid boilerplate property setup.

---

## Extensions

This design can be extended to:

- SEIR models with an `EXPOSED` state and incubation timer
- Multiple patches (nodes) by using a grid `M x N > 1`
- Age-structured populations by adding `date_of_birth` and age-bin logic

Let me know if you want this rewritten as a reusable module or notebook.
