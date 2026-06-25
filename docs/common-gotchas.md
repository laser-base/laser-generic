# Common gotchas

A quick reference to the API details and pitfalls that most often trip people up when
building models with `laser-generic`. Each item shows the correct pattern and the mistake
to avoid (including the error you'll see when you get it wrong).

## Imports and namespace

`laser-generic` is a namespace package — always import it with dot notation:

```python
from laser.generic import Model
from laser.generic.SIR import Susceptible, Infectious, Recovered, Transmission
```

Do **not** use an underscore (`import laser_generic`) — that is not how the package is
imported.

## Disease components live in per-model modules

Import the disease components from the module that matches the model you are building —
`laser.generic.SI`, `.SIS`, `.SIR`, `.SIRS`, `.SEIR`, `.SEIRS`. Each exposes the components
appropriate to that model (`Susceptible`, `Infectious`, `Transmission`, and where relevant
`Recovered` / `Exposed`):

```python
from laser.generic.SIR import Susceptible, Infectious, Recovered, Transmission
```

Prefer these per-model classes. The low-level variants in `laser.generic.components`
(`TransmissionSI`, `TransmissionSE`, `InfectiousIR`, …) have misleading names — e.g.
`components.TransmissionSI` still **requires** an `infdurdist` argument — so reach for them
only if you specifically need them.

### Component constructor arguments differ by model

Build the duration distributions as callables from `laser.core.distributions`
(`exponential`, `poisson`, `constant_int`, `gamma`, …) and pass them positionally. The
signatures differ by model:

| Model | Susceptible | Exposed | Infectious | Recovered | Transmission |
|-------|-------------|---------|------------|-----------|--------------|
| SI    | `(model)` | — | `(model)` | — | `(model)` |
| SIS   | `(model)` | — | `(model, infdurdist)` | — | `(model, infdurdist)` |
| SIR   | `(model)` | — | `(model, infdurdist)` | `(model)` | `(model, infdurdist)` |
| SIRS  | `(model)` | — | `(model, infdurdist, wandurdist)` | `(model, wandurdist)` | `(model, infdurdist)` |
| SEIR  | `(model)` | `(model, expdurdist, infdurdist)` | `(model, infdurdist)` | `(model)` | `(model, expdurdist)` |
| SEIRS | `(model)` | `(model, expdurdist, infdurdist)` | `(model, infdurdist, wandurdist)` | `(model, wandurdist)` | `(model, expdurdist)` |

Key points:

- `Transmission` takes `expdurdist` for the `SE*` models, but `infdurdist` for the rest.
- The `*RS` (waning-immunity) models add a waning-duration distribution `wandurdist`.
- `Susceptible(model)` takes only the model.

```python
from laser.core.distributions import exponential
infdurdist = exponential(scale=8.0)   # mean infectious duration ~8 days
model.components = [Susceptible(model), Infectious(model, infdurdist),
                    Recovered(model), Transmission(model, infdurdist)]
```

## Parameters (`PropertySet`)

The model's `PropertySet` must include `beta` (the per-day transmission rate) and `nticks`
(number of ticks to run). The `Transmission` component reads `model.params.beta` and raises
`AttributeError: 'PropertySet' object has no attribute 'beta'` if it is missing.

The random seed is **optional**: `Model.__init__` looks for `prng_seed`, `prngseed`, or
`seed` (in that order) and falls back to a fixed default of `20260101` for reproducibility.
Set one only when you want to control or vary the seed yourself.

```python
from laser.generic.utils import PropertySet
params = PropertySet({"nticks": 180, "beta": 0.3, "seed": 42})  # seed optional
```

Also include `inf_mean` (mean infectious duration, in days) whenever you call
`seed_infections_randomly` / `seed_infections_in_patch` — they read it to set the infection
timer of seeded agents.

`PropertySet` supports both attribute access (`model.params.beta`) and dict-style subscript
access (`model.params["beta"]`), but **attribute access is strongly preferred** — it reads
more naturally and is the style used throughout the laser-generic codebase. For optional
parameters, use `getattr(model.params, "name", default)`.

## Seeding and importation

`seed_infections_randomly(model, n)` (from `laser.generic.importation`) infects `n` random
susceptibles — it sets their `state` and draws their `itimer` from `model.params.inf_mean`
(so `inf_mean` must be in `params`). Use it for **initial** seeding, **before** `model.run()`.

For **periodic / mid-run importation** it is not sufficient on its own: it does not update the
node-level compartment counts, so `Transmission` won't see the imported infections and the
outbreak never starts (peak stays 0). Instead write a small custom component whose
`step(self, tick)` (a) picks random susceptible agents, (b) sets their `state` to
`INFECTIOUS` and their `itimer` from the infectious-duration distribution, and (c) updates
the next tick's node counts:

```python
self.model.nodes.S[tick + 1] -= n_by_node
self.model.nodes.I[tick + 1] += n_by_node
```

Note: `Infect_Random_Agents` is **not** a per-tick component (it has no `step` method), and
there is no class named `Importation` to import.

## Vital dynamics and age pyramids

`MortalityByCDR(model, mortalityrates)` (deaths) and
`ConstantPopVitalDynamics(model, recycle_rates)` (balanced births and deaths, total
population conserved) need **no** age pyramid.

`BirthsByCBR(model, birthrates, pyramid=...)` **requires** a `pyramid` argument, and that
argument must be an `AliasedDistribution` (not a raw array). `AliasedDistribution` and
`load_pyramid_csv` are imported from **`laser.core.demographics`** (not
`laser.core.distributions`). A simple synthetic pyramid is fine:

```python
import numpy as np
from laser.core.demographics import AliasedDistribution

# counts per yearly age bin (a stable, exponentially-declining age structure)
pyramid = AliasedDistribution(np.array(1000 * np.exp(-0.02 * np.arange(89))))
model.components.append(BirthsByCBR(model, birthrate_map, pyramid=pyramid))
```

Passing a degenerate or empty array as the pyramid raises `ValueError: high <= 0`.

## Reading results

Compartment time series are tracked at the **node level** as arrays on `model.nodes`:
`model.nodes.S`, `model.nodes.I`, `model.nodes.E`, `model.nodes.R` (plus
`newly_infectious` / `newly_infected` where relevant). Each is allocated with length
`nticks + 1`, so the shape is `[nticks + 1, n_nodes]` — there is one row per tick from `0`
through `nticks` inclusive. (This means `I[-1]` is the state at tick `nticks`, not
`nticks - 1`.) There is no `StateTracker` or `ResultsWriter` class.

```python
I = model.nodes.I[:, 0]      # single-patch infectious time series
peak = I.max()
final = I[-1]
total_per_tick = model.nodes.I.sum(axis=1)   # across all patches
```

Reading `model.nodes.I[:, 0]` **after** `model.run()` is always safe. Do not compute
compartment counts by testing the per-agent `model.people.state` against a bare integer
(e.g. `people.state == 1`); if you need per-agent state, compare against the `State` enum
value (e.g. `State.INFECTIOUS.value`, with `from laser.generic import State`).

When **recording a series from inside a custom component's `step(self, tick)`**, the node
array for the *current* tick may not be written yet (the disease components write it during
their own step), so `model.nodes.I[tick, 0]` can read 0. Count agents directly instead:

```python
from laser.generic import State
infectious_now = int((self.model.people.state == State.INFECTIOUS.value).sum())
```

## Writing a custom component

A component is a plain class the model calls as `component.step(tick)` each tick. Store the
model in `__init__` as `self.model = model` and reference `self.model.people` /
`self.model.nodes` inside `step` — never a bare module-level `model`. `model.nodes` is a
`LaserFrame`; the number of patches is `model.nodes.count` (an int) — don't treat
`model.nodes` itself as a number.

```python
class MyComponent:
    def __init__(self, model):
        self.model = model

    def step(self, tick):
        n_patches = self.model.nodes.count
        ...
```

## Spatial models and migration

Build the migration matrix with `gravity` from `laser.generic.model`. Its signature is
`gravity(pops, distances, k, a, b, c)` — **one** 1-D population array and a **2-D** distance
matrix (not two population arrays); then row-normalize:

```python
from laser.generic.model import gravity, distance
m = gravity(pops, distances, k=1.0, a=1.0, b=1.0, c=2.0)
m /= m.sum(axis=1, keepdims=True)
```

Build the N×N inter-patch distance matrix with `distance(lats, lons, lats, lons)` (from
`laser.generic.model`), where `lats` and `lons` are 1-D arrays of patch latitudes and
longitudes. It takes **two coordinate pairs** (source then destination) — pass the same
`lats`/`lons` twice to get the full pairwise great-circle matrix. Calling it with just two
arguments raises a `TypeError`. Do **not** pass a single 2-D coordinate array or call
`distance(centroids)`.

For a synthetic set of patches, define `lats` / `lons` as plain numpy arrays yourself (or read
existing scenario columns with `scenario["lat"].to_numpy()`). Only use `get_centroids` when
you have a real geopandas `GeoDataFrame` of polygons — it returns a geopandas `GeoSeries` of
point centroids (read coordinates via `centroids.x` / `centroids.y`), **not** an N×2 numpy
array or a `(lats, lons)` tuple. Calling `get_centroids(scenario)` on an ordinary pandas
DataFrame raises `AttributeError: 'DataFrame' object has no attribute 'to_crs'`.

## Working with laser-core primitives

`laser-generic` is built on `laser-core`, so two `laser-core` conventions apply throughout:

- **`LaserFrame` properties are modified in place.** Use `frame.S[:] = new_values`; never
  rebind a property (`frame.S = ...` raises `RuntimeError: Cannot reassign property`). Access
  properties by attribute (`frame.state`), not by subscript (`frame['state']`).
- **`PropertySet` prefers attribute access** (`model.params.beta`) over dict-style subscript
  access (`model.params["beta"]`) — both work, but attribute access is the idiomatic style
  throughout the codebase (see *Parameters* above).
