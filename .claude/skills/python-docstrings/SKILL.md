---
name: python-docstrings
description: Guidelines for writing Python docstrings in this project. Always use when writing, updating, editing, or reviewing docstrings for laser-generic components, functions, or classes. Covers Google-style format, MkDocs cross-references, usage examples, and researcher workflow context.
---

# Python docstring guidelines

This project uses **Google-style docstrings** rendered by [mkdocstrings](https://mkdocstrings.github.io/) with the Material for MkDocs theme. Every public class, method, and function should have a docstring that helps a disease-modeling researcher understand not just *what* the object does, but *where it fits* in their workflow.

---

## Structure of a good docstring

```python
def some_function(arg1, arg2):
    """One-line summary ending with a period.

    Two to four sentences of prose explaining what this does and why a researcher
    would reach for it. Mention the disease model context — is this used during
    initialization, per-tick stepping, or post-simulation analysis?

    Link to related objects using MkDocs cross-reference syntax, e.g.
    [`Model`][laser.generic.model.Model] or [`Susceptible`][laser.generic.components.Susceptible].

    Args:
        arg1 (np.ndarray): Description including shape, dtype, and units where
            relevant. E.g. "Shape ``(nticks+1, num_nodes)``, dtype float32."
        arg2 (int): Description. State default values explicitly if they matter.

    Returns:
        np.ndarray: What is returned, its shape, dtype, and interpretation.
            For None returns, omit this section entirely.

    Raises:
        ValueError: When and why this is raised.

    Example:
        Show the most common researcher workflow. Use a self-contained snippet
        that a researcher can run or adapt directly:

        ```python
        model = Model(scenario, params)
        result = some_function(model.nodes.S[0], model.params.nticks)
        ```

    Note:
        Reserve for non-obvious caveats — performance warnings, thread safety,
        or behaviours that diverge from what the name implies.
    """
```

---

## Cross-references to other Python objects

The project uses the [`autorefs`](https://mkdocstrings.github.io/autorefs/) plugin so you can link to any documented Python object directly from docstring prose.

**Syntax:**

| What you write | What renders |
|---|---|
| `` [`Model`][laser.generic.model.Model] `` | clickable link labelled `Model` |
| `` [`Susceptible`][laser.generic.components.Susceptible] `` | link to the Susceptible class |
| `` [`State`][laser.generic.shared.State] `` | link to the State enum |
| `` [`ValuesMap`][laser.generic.utils.ValuesMap] `` | link to ValuesMap |

Use these links in the **description prose** to guide a reader toward related components they will likely need next. Do not stuff links into Args or Returns sections — prose is the right place.

**Common objects to link to:**

- [`Model`][laser.generic.model.Model] — the top-level simulation container
- [`Susceptible`][laser.generic.components.Susceptible] — must appear before transmission components
- [`Exposed`][laser.generic.components.Exposed] — incubation period component for SEIR models
- [`State`][laser.generic.shared.State] — enum of agent health states (SUSCEPTIBLE, EXPOSED, INFECTIOUS, RECOVERED)
- [`ValuesMap`][laser.generic.utils.ValuesMap] — time-varying per-node parameter container

---

## Researcher workflow context

A researcher using this library follows a predictable workflow. Your docstring should anchor the object to one of these stages:

```
1. Prepare scenario     →  GeoDataFrame with population, geometry, initial S/I/E/R counts
2. Set parameters       →  PropertySet with nticks, beta, capacity_safety_factor, prng_seed, ...
3. Build component list →  [Susceptible(model), Transmission(model, ...), ...]
4. Initialise model     →  model = Model(scenario, params)
5. Run simulation       →  model.run()
6. Analyse results      →  model.nodes.S, model.nodes.I, component.plot()
```

In the description prose, state clearly:
- **Which stage** this object belongs to (initialisation, per-tick step, post-run analysis)
- **What comes before it** in the component list, if ordering matters
- **What it produces** that downstream components or analysis code will consume

Example phrasing:
> "Add this component to `model.components` *before* any transmission component. It initialises the `S` array consumed by [`SIR.Transmission`][laser.generic.SIR.Transmission] at every tick."

---

## Usage examples

Every public class and every non-trivial function needs an `Example:` section. A good example:

1. **Shows the full calling context**, not just the function call in isolation.
2. **Uses realistic variable names** (`scenario`, `params`, `model`) that match the tutorials.
3. **Demonstrates the most common researcher use case** — not edge cases.
4. **Is short** — three to eight lines. If you need more, split into two examples.

**Class example (component):**

```python
Example:
    Assemble a basic SIR model by listing components in execution order:

    ```python
    model.components = [
        Susceptible(model),
        SIR.Transmission(model, beta=0.3),
        SIR.Infectious(model, infdurdist, infdurmin=1),
        SIR.Recovered(model),
    ]
    model.run()
    model.nodes.S  # shape (nticks+1, num_nodes)
    ```
```

**Function example:**

```python
Example:
    Compute gravity-weighted migration rates between patches:

    ```python
    rates = gravity(populations, distances, k=1.0, a=1.0, b=2.0)
    migration_matrix = row_normalizer(rates)
    ```
```

---

## Class docstrings

For classes, the docstring lives on the class body (not `__init__`), because `merge_init_into_class: true` is set in `mkdocs.yml`. Structure it as:

```
One-line summary.

Two to four sentences of context: what disease model states this manages,
which stage of the researcher workflow it belongs to, and what it produces.
Link to closely related components.

Args:
    model (Model): The simulation model this component is attached to.
    some_dist (Callable): Brief description including expected signature.
    some_min (int): Minimum value in days. Defaults to 1.
    validating (bool): If True, runs pre/post validation hooks each tick.
        Enable during development; disable for production runs.

Example:
    ...

Note:
    ...
```

Do not repeat the class name in the one-line summary. Write "Manages the
susceptible compartment" not "Susceptible manages the susceptible compartment."

---

## What to omit

- **Do not document private methods** (`_check_flow_vs_census`, `nb_timer_update`, etc.) unless there is a strong reason.
- **Do not add a docstring to `__init__`** — put everything on the class.
- **Do not restate the type annotation** in the description ("arg1 is a numpy array" when the type is already `np.ndarray`). Use the prose to explain *meaning*, not *type*.
- **Do not write "This function/class ..."** — just say what it does directly.
- **Do not add a `Returns:` section for `None`**.
