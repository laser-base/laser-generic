"""
This module provides utility functions for the laser-measles project.

"""

import time
from math import ceil
from typing import Any
from typing import ClassVar

import geopandas as gpd
import numpy as np

from laser.core import PropertySet
from laser.generic.shared import State
from pyproj import Transformer
from shapely.geometry import Point

__all__ = ["TimingStats", "ValuesMap", "get_centroids"]


class ValuesMap:
    """
    A class to efficiently represent values mapped over nodes and time steps.

    Arguments:
        nnodes (int): Number of nodes.
        nticks (int): Number of time steps.

    Methods to create ValuesMap from different data sources:
        - from_scalar(scalar: float, nticks: int, nnodes: int)
        - from_timeseries(data: np.ndarray, nnodes: int)
        - from_nodes(data: np.ndarray, nticks: int)
        - from_array(data: np.ndarray, writeable: bool = False)
    """

    def __init__(self, nnodes: int, nticks: int):
        self._nnodes = nnodes
        self._nticks = nticks

        return

    @staticmethod
    def from_scalar(scalar: float, nticks: int, nnodes: int) -> "ValuesMap":
        """
        Create a ValuesMap with the same scalar value for all nodes and time steps.

        Args:
            scalar (float): The scalar value to fill the map.
            nnodes (int): Number of nodes.
            nticks (int): Number of time steps.

        Returns:
            ValuesMap: The created ValuesMap instance.
        """
        assert scalar >= 0.0, "scalar must be non-negative"
        assert nnodes > 0, "nnodes must be greater than 0"
        assert nticks > 0, "nticks must be greater than 0"
        instance = ValuesMap(nnodes=nnodes, nticks=nticks)
        tmp = np.array([[scalar]], dtype=np.float32)
        instance._data = np.broadcast_to(tmp, (nticks, nnodes))

        return instance

    @staticmethod
    def from_timeseries(data: np.ndarray, nnodes: int, nticks: int = None) -> "ValuesMap":
        """
        Create a ValuesMap from a time series array for all nodes.

        All nodes have the same time series data.

        nticks is inferred from the length of data if not explicitly provided.

        Args:
            data (np.ndarray): 1D array of time series data.
            nnodes (int): Number of nodes.
            nticks (int, optional): Number of ticks. Defaults to None.

        Returns:
            ValuesMap: The created ValuesMap instance.
        """
        assert np.all(data >= 0.0), "data must be non-negative"
        assert len(data.shape) == 1, "data must be a 1D array"
        assert data.shape[0] > 0, "data must have at least one element"
        assert nnodes > 0, "nnodes must be greater than 0"
        nticks = data.shape[0] if nticks is None else nticks
        assert nticks > 0, "nticks must be greater than 0"
        instance = ValuesMap(nnodes=nnodes, nticks=nticks)
        if data.shape[0] == nticks:
            row = data
        else:
            row = data[np.arange(nticks) % data.shape[0]]
        instance._data = np.broadcast_to(row[:, None], (nticks, nnodes))

        return instance

    @staticmethod
    def from_nodes(data: np.ndarray, nticks: int) -> "ValuesMap":
        """
        Create a ValuesMap from a nodes array for all time steps.

        All time steps have the same node data.

        nnodes is inferred from the length of data.

        Args:
            data (np.ndarray): 1D array of node data.
            nticks (int): Number of time steps.

        Returns:
            ValuesMap: The created ValuesMap instance.
        """
        assert np.all(data >= 0.0), "data must be non-negative"
        assert len(data.shape) == 1, "data must be a 1D array"
        assert data.shape[0] > 0, "data must have at least one element"
        assert nticks > 0, "nticks must be greater than 0"
        nnodes = data.shape[0]
        instance = ValuesMap(nnodes=nnodes, nticks=nticks)
        instance._data = np.broadcast_to(data[None, :], (nticks, nnodes))

        return instance

    @staticmethod
    def from_array(data: np.ndarray, writeable: bool = False) -> "ValuesMap":
        """
        Create a ValuesMap from a 2D array of data.

        Args:
            data (np.ndarray): 2D array of shape (nticks, nnodes).
            writeable (bool): If True, the underlying data array is writeable and can be modified during simulation. Default is False.

        Returns:
            ValuesMap: The created ValuesMap instance.
        """
        assert np.all(data >= 0.0), "data must be non-negative"
        assert len(data.shape) == 2, "data must be a 2D array"
        assert data.shape[0] > 0, "data must have at least one row"
        assert data.shape[1] > 0, "data must have at least one column"
        nticks, nnodes = data.shape
        instance = ValuesMap(nnodes=nnodes, nticks=nticks)
        instance._data = data.astype(np.float32)
        instance._data.flags.writeable = writeable

        return instance

    @property
    def nnodes(self):
        """Number of nodes."""
        return self._nnodes

    @property
    def nticks(self):
        """Number of time steps."""
        return self._nticks

    @property
    def shape(self):
        """Shape of the underlying data array (nticks, nnodes)."""
        return self._data.shape

    @property
    def values(self):
        """Underlying data array of shape (nticks, nnodes)."""
        return self._data

    def __getitem__(self, access):
        return self._data[access]


class TimingContext:
    """Internal class for timing context management."""

    def __init__(self, label: str, stats: "TimingStats", parent: dict) -> None:  # type: ignore
        self.label = label
        self.stats = stats
        self.parent = parent
        self.children = {}
        self.ncalls = 0
        self.elapsed = 0
        self.start = 0
        self.end = 0

        return

    def __enter__(self):
        self.ncalls += 1
        self.stats._enter(self)
        self.start = time.perf_counter_ns()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.perf_counter_ns()
        self.elapsed += self.end - self.start
        self.stats._exit(self)

        return

    @property
    def inclusive(self) -> int:
        return self.elapsed

    @property
    def exclusive(self) -> int:
        excl = self.elapsed
        for child in self.children.values():
            excl -= child.elapsed

        return excl


class _TimingStats:
    """
    Internal class for managing timing statistics.
    """

    def __init__(self) -> None:
        self.frozen = False
        self.context = {}
        self.root = self.start("root")
        self.root.__enter__()

        return

    def start(self, label: str) -> TimingContext:
        """Create a timing context with the given label."""
        assert self.frozen is False

        if label not in self.context:
            self.context[label] = TimingContext(label, self, self.context)

        return self.context[label]

    def _enter(self, context: TimingContext) -> None:
        self.context = context.children
        return

    def _exit(self, context: TimingContext) -> None:
        assert self.context is context.children
        self.context = context.parent
        return

    def freeze(self) -> None:
        """Freeze the timing statistics."""
        assert self.frozen is False
        self.root.__exit__(None, None, None)
        self.frozen = True

        return

    _scale_factors: ClassVar[dict[str, float]] = {
        "ns": 1,
        "nanoseconds": 1,
        "us": 1e3,
        "µs": 1e3,
        "microseconds": 1e3,
        "ms": 1e6,
        "milliseconds": 1e6,
        "s": 1e9,
        "sec": 1e9,
        "seconds": 1e9,
    }

    def to_string(self, scale: str = "ms") -> str:
        assert self.frozen is True

        assert scale in self._scale_factors
        factor = self._scale_factors[scale]

        lines = []

        def _recurse(node: TimingContext, depth: int) -> None:
            indent = "    " * depth
            tot_time = node.elapsed / factor
            avg_time = node.elapsed / node.ncalls / factor if node.ncalls > 0 else 0
            exc_time = node.exclusive / factor
            lines.append(
                f"{indent}{node.label}: {node.ncalls} calls, total {tot_time:.3f} {scale}, avg {avg_time:.3f} {scale}, excl {exc_time:.3f} {scale}"
            )
            for child in node.children.values():
                _recurse(child, depth + 1)

            return

        _recurse(self.root, 0)
        return "\n".join(lines)

    def to_dict(self, scale: str = "ms") -> dict:
        assert self.frozen is True

        assert scale in self._scale_factors
        factor = self._scale_factors[scale]

        def _recurse(node: TimingContext) -> dict:
            result = {
                "label": node.label,
                "ncalls": node.ncalls,
                "inclusive_ns": node.inclusive / factor,
                # "exclusive_ns": node.exclusive / factor,
                "children": [],
            }
            for child in node.children.values():
                result["children"].append(_recurse(child))

            return result

        return _recurse(self.root)


TimingStats = _TimingStats()


def validate(pre, post):
    """
    Decorator to add pre- and post-validation to a method.

    Calls the given pre- and post-validation methods if the model or component is in validating mode.
    """

    def decorator(func):
        def wrapper(self, tick: int, *args, **kwargs):
            if pre and (getattr(self.model, "validating", False) or getattr(self, "validating", False)):
                with TimingStats.start(pre.__name__):
                    getattr(self, pre.__name__)(tick)
            result = func(self, tick, *args, **kwargs)
            if post and (getattr(self.model, "validating", False) or getattr(self, "validating", False)):
                with TimingStats.start(post.__name__):
                    getattr(self, post.__name__)(tick)
            return result

        return wrapper

    return decorator


def get_centroids(gdf: gpd.GeoDataFrame) -> np.ndarray:
    """Get centroids of geometries in gdf in degrees (EPSG:4326)."""

    gdf_3857 = gdf.to_crs(epsg=3857)
    centroids_3857 = gdf_3857.geometry.centroid

    # centroids_3857.to_crs(epsg=4326) emits a Warning is there is only one point (one node)
    if len(centroids_3857) > 1:
        centroids_deg = centroids_3857.to_crs(epsg=4326)
    else:
        # Explicitly transform the single centroid
        transformer = Transformer.from_crs(3857, 4326, always_xy=True)
        x, y = centroids_3857.x.values[0], centroids_3857.y.values[0]
        lon, lat = transformer.transform(x, y)
        centroids_deg = gpd.GeoSeries([Point(lon, lat)], crs="EPSG:4326")

    return centroids_deg


# Want to think about the ways to seed infections.  Not all infections have a timer!
def seed_infections_randomly(model: Any, ninfections: int = 100) -> np.ndarray:
    """
    Randomly seed initial infections across the entire people.

    This function selects up to `ninfections` susceptible individuals at random
    from the full people. It marks them as infected by:
    - Setting their infection timer (`itimer`) to the model's mean infectious duration (`inf_mean`),
    - Setting their susceptibility to zero.

    Args:
        model: The simulation model, which must contain a `people` with
               `susceptibility`, `itimer`, and `nodeid` arrays, and a `params` object with `inf_mean`.
        ninfections (int, optional): The number of individuals to infect. Defaults to 100.

    Returns:
        np.ndarray: The node IDs of the newly infected individuals.
    """

    # Seed initial infections in random locations at the start of the simulation
    pop = model.people
    params = model.params

    myinds = np.flatnonzero(pop.state == State.SUSCEPTIBLE.value)
    if len(myinds) > ninfections:
        myinds = np.random.permutation(myinds)[:ninfections]

    pop.itimer[myinds] = params.inf_mean
    pop.state[myinds] = State.INFECTIOUS.value
    inf_nodeids = pop.nodeid[myinds]

    return inf_nodeids


def seed_infections_in_patch(model: Any, ipatch: int, ninfections: int = 1) -> None:
    """
    Seed initial infections in a specific patch of the people at the start of the simulation.

    This function randomly selects up to `ninfections` individuals from the specified patch
    who are currently susceptible (state == State.SUSCEPTIBLE) and marks them as infected by:
      - Setting their infection timer (`itimer`) to the model's mean infectious duration (`inf_mean`),
      - Setting their infection `state` to State.INFECTIOUS.

    Args:
        model: The simulation model containing the people and parameters. It must expose:
               - model.people.state (integer infection state),
               - model.people.itimer (infection timers),
               - model.people.nodeid (patch index),
               - model.params.inf_mean (mean infectious period).
        ipatch (int): The identifier of the patch where infections should be seeded.
        ninfections (int, optional): The number of initial infections to seed. Defaults to 1.

    Returns:
        None
    """
    pop = model.people

    # Candidates: susceptible individuals in the target patch
    susceptible_in_patch = np.where((pop.state == State.SUSCEPTIBLE.value) & (pop.nodeid == ipatch))[0]

    if len(susceptible_in_patch) == 0:
        # Nothing to do: no susceptibles in this patch
        return

    # If there are more candidates than requested infections, sample without replacement
    if len(susceptible_in_patch) > ninfections:
        susceptible_in_patch = np.random.choice(susceptible_in_patch, ninfections, replace=False)

    # Set timers and mark as infectious
    pop.itimer[susceptible_in_patch] = model.params.inf_mean
    pop.state[susceptible_in_patch] = State.INFECTIOUS.value

    return


def get_default_parameters() -> PropertySet:
    """
    Returns a default PropertySet with common parameters used across laser-generic models.

    Each parameter in the returned PropertySet is described below, along with its default value and rationale:

        nticks (int, default=730): Number of simulation ticks (days). Default is 2 years (365*2), which is a typical duration for seasonal epidemic simulations.
        beta (float, default=0.15): Transmission rate per contact. Chosen as a moderate value for SIR-type models to reflect realistic disease spread.
        biweekly_beta_scalar (list of float, default=[1.0]*biweekly_steps): Scalar for beta for each biweekly period. Default is 1.0 for all periods, meaning no seasonal variation unless specified.
        cbr (float, default=0.03): Constant birth rate. Set to 0.03 to represent a typical annual birth rate in people models.
        exp_shape (float, default=2.0): Shape parameter for the exposed period distribution. Default chosen for moderate dispersion.
        exp_scale (float, default=2.0): Scale parameter for the exposed period distribution. Default chosen for moderate mean duration.
        inf_mean (float, default=4.0): Mean infectious period (days). Set to 4.0 to reflect typical infectious durations for diseases like measles.
        inf_sigma (float, default=1.0): Standard deviation of infectious period. Default is 1.0 for moderate variability.
        seasonality_factor (float, default=0.2): Amplitude of seasonal forcing. Chosen to allow moderate seasonal variation in transmission.
        seasonality_phase (float, default=0.0): Phase offset for seasonality. Default is 0.0, meaning no phase shift.
        importation_count (int, default=1): Number of cases imported per importation event. Default is 1 for sporadic importation.
        importation_period (int, default=30): Days between importation events. Default is 30 to represent monthly importation.
        importation_start (int, default=0): Start day for importation events. Default is 0 (simulation start).
        importation_end (int, default=730): End day for importation events. Default is 2 years (365*2).
        seed (int, default=123): Random seed for reproducibility. Default is 123.
        verbose (bool, default=False): If True, enables verbose output. Default is False for minimal output.
    These values are chosen to be broadly reasonable for seasonal SIR-type models with importation.

    We need a function like this because even-though laser-core requires no particular param name,
    laser-generic code does presume certain parameters and there's no elegant way to just discover
    what those are. So we put them here.
    """
    nticks = 365 * 2
    biweekly_steps = ceil(nticks / 14)
    return PropertySet(
        {
            "nticks": nticks,
            "beta": 0.15,
            "biweekly_beta_scalar": [1.0] * biweekly_steps,
            "cbr": 0.03,
            "exp_shape": 2.0,
            "exp_scale": 2.0,
            "inf_mean": 4.0,
            "inf_sigma": 1.0,
            "seasonality_factor": 0.2,
            "seasonality_phase": 0.0,
            "importation_count": 1,
            "importation_period": 30,
            "importation_start": 0,
            "importation_end": 365 * 2,
            "prng_seed": 314159265,
            "verbose": False,
        }
    )
