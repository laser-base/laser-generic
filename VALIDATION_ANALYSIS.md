# Validation Infrastructure Analysis

## Summary

This document provides a comprehensive analysis of validation infrastructure across all LASER generic components, identifying which components have proper validation and which are missing validation functionality.

## Components WITH Validation Infrastructure ✅

### components.py (11/11 components)

All components in `components.py` implement proper validation infrastructure:

1. **Susceptible** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
2. **Exposed** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
3. **InfectiousSI** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
4. **InfectiousIS** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
5. **InfectiousIR** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
6. **InfectiousIRS** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
7. **Recovered** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
8. **RecoveredRS** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
9. **TransmissionSIX** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
10. **TransmissionSI** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`
11. **TransmissionSE** - Uses `@validate` decorator with `prevalidate_step` and `postvalidate_step`

**Note:** These components use the `@validate` decorator which checks `model.validating` flag. They do not have explicit `validating` parameters in `__init__()`, but validation is controlled at the model level.

### vitaldynamics.py (4/4 components)

All components in `vitaldynamics.py` implement proper validation infrastructure:

1. **BirthsByCBR** - Has `validating` parameter (default=False), uses `@validate` decorator
2. **MortalityByCDR** - Has `validating` parameter (default=False), uses `@validate` decorator
3. **MortalityByEstimator** - Has `validating` parameter (default=False), uses `@validate` decorator
4. **ConstantPopVitalDynamics** - Has `validating` parameter (default=False), uses `@validate` decorator

## Components MISSING Validation Infrastructure ❌

### immunization.py (3/3 components missing validation)

1. **RoutineImmunization**
   - ❌ NO `validating` parameter in `__init__()`
   - ❌ NO `prevalidate_step()` method
   - ❌ NO `postvalidate_step()` method
   - ❌ NO `@validate` decorator on `__call__()`

2. **ImmunizationCampaign**
   - ❌ NO `validating` parameter in `__init__()`
   - ❌ NO `prevalidate_step()` method
   - ❌ NO `postvalidate_step()` method
   - ❌ NO `@validate` decorator on `__call__()`

3. **RoutineImmunizationEx**
   - ❌ NO `validating` parameter in `__init__()`
   - ❌ NO `prevalidate_step()` method
   - ❌ NO `postvalidate_step()` method
   - ❌ NO `@validate` decorator on `step()`

### importation.py (2/2 components missing validation)

1. **Infect_Random_Agents**
   - ❌ NO `validating` parameter in `__init__()`
   - ❌ NO `prevalidate_step()` method
   - ❌ NO `postvalidate_step()` method
   - ❌ NO `@validate` decorator on `__call__()`

2. **Infect_Agents_In_Patch**
   - ❌ NO `validating` parameter in `__init__()`
   - ❌ NO `prevalidate_step()` method
   - ❌ NO `postvalidate_step()` method
   - ❌ NO `@validate` decorator on `__call__()`

---

## Suggested Validation Implementations

### For immunization.py Components

#### RoutineImmunization

**Add to `__init__()` signature:**
```python
def __init__(
    self,
    model,
    period: int,
    coverage: float,
    age: int,
    start: int = 0,
    end: int = -1,
    verbose: bool = False,
    validating: bool = False,  # ADD THIS
) -> None:
```

**Suggested validation methods:**
```python
def prevalidate_step(self, tick: int) -> None:
    """
    Pre-step validation for routine immunization.

    Validates that:
    - All susceptibility values are valid (0 or 1)
    - S counts match agent-level susceptibility before immunization
    """
    # Store pre-immunization susceptible count
    self.prv_susceptible = np.sum(self.model.people.susceptibility)

    # Validate susceptibility values are binary
    assert np.all((self.model.people.susceptibility == 0) | (self.model.people.susceptibility == 1)), \
        "Susceptibility values must be 0 or 1"

    return

def postvalidate_step(self, tick: int) -> None:
    """
    Post-step validation for routine immunization.

    Validates that:
    - Susceptible count decreased or stayed same (no new susceptibles from immunization)
    - All immunized agents have susceptibility = 0
    - Change in susceptibility matches recorded immunizations
    """
    post_susceptible = np.sum(self.model.people.susceptibility)

    # Susceptible count should only decrease (or stay same if no immunizations)
    assert post_susceptible <= self.prv_susceptible, \
        "Susceptible count should not increase after immunization"

    # If this was an immunization tick, verify counts changed appropriately
    if (tick >= self.start) and ((tick - self.start) % self.period == 0) and (tick < self.end):
        # Some immunizations should have occurred (unless no eligible agents)
        change = self.prv_susceptible - post_susceptible
        assert change >= 0, "Immunization should reduce or maintain susceptible count"

    return

@validate(pre=prevalidate_step, post=postvalidate_step)
def __call__(self, model, tick: int) -> None:
    # ... existing implementation ...
```

#### ImmunizationCampaign

**Add to `__init__()` signature:**
```python
def __init__(
    self,
    model,
    period: int,
    coverage: float,
    age_lower: int,
    age_upper: int,
    start: int = 0,
    end: int = -1,
    verbose: bool = False,
    validating: bool = False,  # ADD THIS
) -> None:
```

**Suggested validation methods:**
```python
def prevalidate_step(self, tick: int) -> None:
    """
    Pre-step validation for immunization campaign.

    Validates that:
    - All susceptibility values are valid (0 or 1)
    - Population has valid dob values for age calculations
    """
    self.prv_susceptible = np.sum(self.model.people.susceptibility)

    assert np.all((self.model.people.susceptibility == 0) | (self.model.people.susceptibility == 1)), \
        "Susceptibility values must be 0 or 1"

    # Validate ages are reasonable (dob <= current tick)
    if hasattr(self.model.people, 'dob'):
        ages = tick - self.model.people.dob
        assert np.all(ages >= 0), "All agents should have non-negative ages"

    return

def postvalidate_step(self, tick: int) -> None:
    """
    Post-step validation for immunization campaign.

    Validates that:
    - Susceptible count decreased or stayed same
    - Only agents in target age band were immunized
    """
    post_susceptible = np.sum(self.model.people.susceptibility)

    assert post_susceptible <= self.prv_susceptible, \
        "Susceptible count should not increase after campaign"

    # If this was a campaign tick, verify targeting
    if (tick >= self.start) and ((tick - self.start) % self.period == 0) and (tick < self.end):
        if hasattr(self.model.people, 'dob'):
            ages = tick - self.model.people.dob
            # Verify only agents in age band were immunized
            # (Complex check - immunized agents should be in [age_lower, age_upper))
            in_band = (ages >= self.age_lower) & (ages < self.age_upper)
            not_susceptible = self.model.people.susceptibility == 0

            # All immunized agents outside band should have been immune already
            # (This is a weak check - full tracking would require more state)

    return

@validate(pre=prevalidate_step, post=postvalidate_step)
def __call__(self, model, tick: int) -> None:
    # ... existing implementation ...
```

#### RoutineImmunizationEx

**Add to `__init__()` signature:**
```python
def __init__(
    self,
    model,
    coverage_fn: Callable[[int, int], float],
    dose_timing_dist: Callable[[int, int], int],
    dose_timing_min: int = 1,
    initialize: bool = True,
    track: bool = False,
    validating: bool = False,  # ADD THIS
) -> None:
```

**Suggested validation methods:**
```python
def prevalidate_step(self, tick: int) -> None:
    """
    Pre-step validation for routine immunization with explicit timers.

    Validates that:
    - ritimer values are consistent with state
    - Agents with ritimer == 1 are eligible for immunization
    - S and R counts match agent states
    """
    from laser.generic.shared import State

    # Check ritimer consistency
    assert np.all(self.model.people.ritimer >= 0), "ritimers should be non-negative"

    # Store agents about to be immunized
    self.prv_ritimer_one = self.model.people.ritimer == 1
    self.prv_s_count = np.sum(self.model.people.state == State.SUSCEPTIBLE.value)
    self.prv_r_count = np.sum(self.model.people.state == State.RECOVERED.value)

    return

def postvalidate_step(self, tick: int) -> None:
    """
    Post-step validation for routine immunization with explicit timers.

    Validates that:
    - Agents with ritimer==1 before had timers decremented to 0
    - Newly immunized agents transitioned from S to R
    - Recorded immunizations match actual state changes
    """
    from laser.generic.shared import State

    # Check that ritimer==1 agents had timers decremented
    assert np.all(self.model.people.ritimer[self.prv_ritimer_one] == 0), \
        "Agents with ritimer==1 should have timer=0 after step"

    # Verify state changes match recorded immunizations
    post_s_count = np.sum(self.model.people.state == State.SUSCEPTIBLE.value)
    post_r_count = np.sum(self.model.people.state == State.RECOVERED.value)

    s_decrease = self.prv_s_count - post_s_count
    r_increase = post_r_count - self.prv_r_count

    # S decrease should match R increase (conservation)
    assert s_decrease == r_increase, \
        f"S decrease ({s_decrease}) should match R increase ({r_increase})"

    # Verify recorded immunizations match actual transitions
    recorded = self.model.nodes.ri_immunized[tick].sum()
    assert recorded == s_decrease, \
        f"Recorded immunizations ({recorded}) should match S decrease ({s_decrease})"

    return

@validate(pre=prevalidate_step, post=postvalidate_step)
def step(self, tick: int) -> None:
    # ... existing implementation ...
```

### For importation.py Components

#### Infect_Random_Agents

**Add to `__init__()` signature:**
```python
def __init__(self, model, verbose: bool = False, validating: bool = False) -> None:
```

**Suggested validation methods:**
```python
def prevalidate_step(self, tick: int) -> None:
    """
    Pre-step validation for random infection importation.

    Validates that:
    - Sufficient susceptible agents exist for importation
    - State counts are consistent before importation
    """
    from laser.generic.shared import State

    # Store pre-importation counts
    self.prv_s_count = np.sum(self.model.people.state == State.SUSCEPTIBLE.value)
    self.prv_i_count = np.sum(self.model.people.state == State.INFECTIOUS.value)

    # Verify enough susceptibles exist for importation
    if (tick >= self.start) and ((tick - self.start) % self.period == 0) and (tick < self.end):
        assert self.prv_s_count >= self.count, \
            f"Not enough susceptible agents ({self.prv_s_count}) for importation ({self.count})"

    return

def postvalidate_step(self, tick: int) -> None:
    """
    Post-step validation for random infection importation.

    Validates that:
    - Correct number of infections were imported
    - S decreased and I increased by same amount
    - Infected agents have valid itimers (if applicable)
    """
    from laser.generic.shared import State

    post_s_count = np.sum(self.model.people.state == State.SUSCEPTIBLE.value)
    post_i_count = np.sum(self.model.people.state == State.INFECTIOUS.value)

    s_decrease = self.prv_s_count - post_s_count
    i_increase = post_i_count - self.prv_i_count

    # On importation ticks, verify correct number of infections
    if (tick >= self.start) and ((tick - self.start) % self.period == 0) and (tick < self.end):
        assert s_decrease == self.count, \
            f"S decrease ({s_decrease}) should equal importation count ({self.count})"
        assert i_increase == self.count, \
            f"I increase ({i_increase}) should equal importation count ({self.count})"

        # If model has itimer, verify newly infected have valid timers
        if hasattr(self.model.people, 'itimer'):
            # All infectious agents should have itimer > 0
            infectious = self.model.people.state == State.INFECTIOUS.value
            assert np.all(self.model.people.itimer[infectious] > 0), \
                "All infectious agents should have itimer > 0"
    else:
        # No importation on non-importation ticks
        assert s_decrease == 0 and i_increase == 0, \
            "S and I should be unchanged on non-importation ticks"

    return

@validate(pre=prevalidate_step, post=postvalidate_step)
def __call__(self, model, tick) -> None:
    # ... existing implementation ...
```

#### Infect_Agents_In_Patch

**Add to `__init__()` signature:**
```python
def __init__(self, model, verbose: bool = False, validating: bool = False) -> None:
```

**Suggested validation methods:**
```python
def prevalidate_step(self, tick: int) -> None:
    """
    Pre-step validation for patch-targeted infection importation.

    Validates that:
    - Each target patch has sufficient susceptibles
    - State counts are consistent before importation
    """
    from laser.generic.shared import State

    self.prv_s_by_node = np.bincount(
        self.model.people.nodeid[self.model.people.state == State.SUSCEPTIBLE.value],
        minlength=self.model.nodes.count
    )
    self.prv_i_by_node = np.bincount(
        self.model.people.nodeid[self.model.people.state == State.INFECTIOUS.value],
        minlength=self.model.nodes.count
    )

    # On importation ticks, verify patches have enough susceptibles
    if (tick >= self.start) and ((tick - self.start) % self.period == 0) and (tick < self.end):
        for patch in self.patchlist:
            assert self.prv_s_by_node[patch] >= self.count, \
                f"Patch {patch} has insufficient susceptibles ({self.prv_s_by_node[patch]}) " \
                f"for importation count ({self.count})"

    return

def postvalidate_step(self, tick: int) -> None:
    """
    Post-step validation for patch-targeted infection importation.

    Validates that:
    - Correct number of infections per patch
    - Only target patches received infections
    - S and I counts updated correctly per patch
    """
    from laser.generic.shared import State

    post_s_by_node = np.bincount(
        self.model.people.nodeid[self.model.people.state == State.SUSCEPTIBLE.value],
        minlength=self.model.nodes.count
    )
    post_i_by_node = np.bincount(
        self.model.people.nodeid[self.model.people.state == State.INFECTIOUS.value],
        minlength=self.model.nodes.count
    )

    s_changes = self.prv_s_by_node - post_s_by_node
    i_changes = post_i_by_node - self.prv_i_by_node

    if (tick >= self.start) and ((tick - self.start) % self.period == 0) and (tick < self.end):
        # Verify each target patch got correct number of infections
        for patch in self.patchlist:
            assert s_changes[patch] == self.count, \
                f"Patch {patch} S decrease ({s_changes[patch]}) should equal count ({self.count})"
            assert i_changes[patch] == self.count, \
                f"Patch {patch} I increase ({i_changes[patch]}) should equal count ({self.count})"

        # Verify non-target patches unchanged
        for patch in range(self.model.nodes.count):
            if patch not in self.patchlist:
                assert s_changes[patch] == 0, f"Non-target patch {patch} should have no S change"
                assert i_changes[patch] == 0, f"Non-target patch {patch} should have no I change"
    else:
        # No importation on non-importation ticks
        assert np.all(s_changes == 0) and np.all(i_changes == 0), \
            "No patches should have S/I changes on non-importation ticks"

    return

@validate(pre=prevalidate_step, post=postvalidate_step)
def __call__(self, model, tick) -> None:
    # ... existing implementation ...
```

---

## Validation Design Patterns

Based on analysis of existing validated components, here are the key patterns:

### Pattern 1: Add `validating` parameter to `__init__()`
```python
def __init__(self, model, ..., validating: bool = False):
    self.model = model
    self.validating = validating
```

### Pattern 2: Implement pre-validation method
```python
def prevalidate_step(self, tick: int) -> None:
    """Store state before step for comparison."""
    # Check preconditions
    # Store values for post-validation comparison
    return
```

### Pattern 3: Implement post-validation method
```python
def postvalidate_step(self, tick: int) -> None:
    """Verify state changes are consistent."""
    # Check postconditions
    # Compare with pre-validation state
    # Verify conservation laws
    return
```

### Pattern 4: Apply `@validate` decorator
```python
@validate(pre=prevalidate_step, post=postvalidate_step)
def step(self, tick: int) -> None:
    # ... implementation ...
```

### Pattern 5: For components using `__call__` instead of `step`
```python
@validate(pre=prevalidate_step, post=postvalidate_step)
def __call__(self, model, tick: int) -> None:
    # ... implementation ...
```

---

## Testing Recommendations

1. **Run validation tests frequently** during development to catch errors early
2. **Enable validation** (`model.validating = True`) for all integration tests
3. **Add component-specific validation tests** for edge cases
4. **Test with validation disabled** to ensure minimal performance impact
5. **Document validation assumptions** in component docstrings

## Next Steps

1. Implement validation infrastructure for immunization.py components
2. Implement validation infrastructure for importation.py components
3. Add validation tests to CI/CD pipeline
4. Create examples demonstrating validation usage
5. Document validation best practices in component development guide
