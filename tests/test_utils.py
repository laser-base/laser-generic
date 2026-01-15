import unittest
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
from laser.generic import utils
from laser.generic.shared import State


class DummyModel:
    """
    Mock model object for testing utility functions without full Model infrastructure.

    WHAT IT PROVIDES:
    - Lightweight test double that mimics the structure of a real Model
    - Provides people, params, and prng attributes needed by utility functions
    - Allows testing seeding functions, validation decorators, and other utilities in isolation

    STRUCTURE:
    - DummyModel.people: Mock people with susceptibility, itimer, nodeid arrays (not used)
    - DummyModel.Params: Mock parameters (inf_mean)
    - DummyModel.PRNG: Mock random number generator with fixed seed (42) for reproducibility

    USAGE:
    Initialize with a people count (default 10), then pass to utility functions
    that expect a model object. Useful for unit testing without running full simulations.
    """

    class people:
        """
        Mock people container for testing utility functions.

        ATTRIBUTES:
        - count: Number of individuals in people
        - susceptibility: Array of susceptibility values (1.0 = susceptible, 0.0 = immune/infected)
        - itimer: Array of infection timers (time remaining in infectious period)
        - nodeid: Array of node assignments (which spatial node each individual belongs to)

        INITIAL STATE:
        All individuals start fully susceptible (susceptibility=1.0) with no infection
        timers (itimer=0.0), distributed across nodes with sequential IDs.
        """

        def __init__(self, count):
            self.count = count
            self.state = np.zeros(count, dtype=np.int8)
            self.itimer = np.zeros(count)
            self.nodeid = np.arange(count)

    class Params:
        """
        Mock parameter container for testing utility functions.

        ATTRIBUTES:
        - inf_mean: Mean infectious period duration (default 5.0 days)

        Used by seeding functions to set initial infection timer values.
        """

        def __init__(self):
            self.inf_mean = 5.0

    class PRNG:
        """
        Mock pseudo-random number generator with deterministic seed for reproducible testing.

        ATTRIBUTES:
        - rng: NumPy random generator initialized with seed=42

        METHODS:
        - integers(low, high): Generate random integers in range [low, high)
        - choice(a, size, replace): Sample from array without/with replacement

        The fixed seed ensures test reproducibility - same test always produces same random values.
        """

        def __init__(self):
            self.rng = np.random.default_rng(42)

        def integers(self, low, high):
            return self.rng.integers(low, high)

        def choice(self, a, size, replace=False):
            return self.rng.choice(a, size, replace=replace)

    def __init__(self, count=10):
        self.people = DummyModel.people(count)
        self.params = DummyModel.Params()
        self.prng = DummyModel.PRNG()
        self.model = self
        self.validating = False


class TestValuesMap(unittest.TestCase):
    """
    Test suite for ValuesMap factory methods that create spatiotemporal parameter arrays.

    ValuesMap is a core utility for representing parameters that vary across time (nticks)
    and space (nnodes). These tests validate the four factory methods for creating ValuesMap
    objects from different input patterns.

    WHAT IS TESTED:
    - from_scalar(): Constant value across all time and space
    - from_timeseries(): Time-varying value, same across all nodes
    - from_nodes(): Spatially-varying value, constant across time
    - from_array(): Fully custom spatiotemporal pattern
    - Replication: Timeseries cycling/repeating when nticks > array length
    """

    def test_from_scalar(self):
        """
        Validate ValuesMap.from_scalar() creates constant values across time and space.

        SCENARIO SETUP:
        - Scalar value: 2.0
        - Shape: 3 nodes × 4 ticks

        WHAT IS TESTED:
        - Output shape matches (nticks, nnodes) = (4, 3)
        - All values in the array equal the input scalar (2.0)

        FAILURE MEANING:
        If this test fails, from_scalar() isn't correctly broadcasting the scalar to fill
        the entire (nticks, nnodes) array. This would break scenarios requiring constant
        parameters like baseline transmission rates or uniform seasonality.
        """
        vm = utils.ValuesMap.from_scalar(2.0, nticks=4, nnodes=3)
        self.assertEqual(vm.shape, (4, 3))
        self.assertTrue(np.all(vm.values == 2.0))

    def test_from_timeseries(self):
        """
        Validate ValuesMap.from_timeseries() replicates time series across all nodes.

        SCENARIO SETUP:
        - Time series: [1, 2, 3] (3 ticks)
        - Nodes: 2
        - Expected shape: (3, 2)

        WHAT IS TESTED:
        - Output shape matches (nticks, nnodes) where nticks = input array length
        - Each timestep has the same value across all nodes (node axis is replicated)
        - Values follow input pattern: tick 0 = 1, tick 1 = 2, tick 2 = 3

        FAILURE MEANING:
        If this test fails, from_timeseries() isn't correctly replicating temporal values
        across the node dimension. This would break temporal seasonality scenarios where
        all locations experience the same time-varying forcing (e.g., summer/winter cycles).
        """
        arr = np.array([1, 2, 3], dtype=np.float32)
        # nticks inferred from array length
        vm = utils.ValuesMap.from_timeseries(arr, nnodes=2)
        self.assertEqual(vm.shape, (3, 2))
        self.assertTrue(np.all(vm.values[0] == 1))
        self.assertTrue(np.all(vm.values[1] == 2))
        self.assertTrue(np.all(vm.values[2] == 3))

    def test_from_nodes(self):
        """
        Validate ValuesMap.from_nodes() replicates spatial pattern across all timesteps.

        SCENARIO SETUP:
        - Node values: [1, 2, 3] (3 nodes)
        - Timesteps: 2
        - Expected shape: (2, 3)

        WHAT IS TESTED:
        - Output shape matches (nticks, nnodes) where nnodes = input array length
        - Each node has the same value across all timesteps (time axis is replicated)
        - Values follow input pattern: node 0 = 1, node 1 = 2, node 2 = 3

        FAILURE MEANING:
        If this test fails, from_nodes() isn't correctly replicating spatial values across
        the time dimension. This would break spatial seasonality scenarios where different
        locations have different but constant forcing (e.g., climate zones with persistent
        differences in transmission favorability).
        """
        arr = np.array([1, 2, 3], dtype=np.float32)
        # nnodes inferred from array length
        vm = utils.ValuesMap.from_nodes(arr, nticks=2)
        self.assertEqual(vm.shape, (2, 3))
        self.assertTrue(np.all(vm.values[0] == arr))
        self.assertTrue(np.all(vm.values[1] == arr))

    def test_from_array(self):
        """
        Validate ValuesMap.from_array() wraps an existing 2D array without modification.

        SCENARIO SETUP:
        - Input array: 2×3 array of ones
        - Shape inferred from array: (2, 3) representing (nticks=2, nnodes=3)

        WHAT IS TESTED:
        - Output shape matches input array shape
        - All values preserved exactly as provided (no replication or transformation)

        FAILURE MEANING:
        If this test fails, from_array() is modifying or incorrectly wrapping the input
        array. This would break fully custom spatiotemporal patterns where users provide
        explicit values for each (time, node) combination (e.g., empirically measured
        transmission rates or complex intervention schedules).
        """
        arr = np.ones((2, 3), dtype=np.float32)
        # nnodes and nticks inferred from array shape
        vm = utils.ValuesMap.from_array(arr)
        self.assertEqual(vm.shape, (2, 3))
        self.assertTrue(np.all(vm.values == 1.0))

    def test_from_timeseries_replicate_full(self):
        """
        Validate ValuesMap.from_timeseries() cycles/repeats pattern when nticks exceeds array length.

        SCENARIO SETUP:
        - Time series: [0, 1, 2, ..., 9] (10-tick pattern)
        - Requested nticks: 50 (5× longer than pattern)
        - Nodes: 2

        WHAT IS TESTED:
        - Output shape is (50, 2) as requested
        - Pattern repeats exactly 5 times: ticks 0-9 use pattern, 10-19 repeat, etc.
        - Modulo arithmetic correctly maps tick i to pattern[i % 10]

        FAILURE MEANING:
        If this test fails, from_timeseries() isn't correctly cycling the pattern to fill
        longer simulations. This would break seasonal forcing scenarios where a yearly pattern
        (e.g., 365 days) needs to repeat across multi-year simulations (e.g., 1095 days/3 years).
        The pattern must tile seamlessly without gaps or discontinuities.
        """
        arr = np.arange(10, dtype=np.float32)
        nticks = 50
        nnodes = 2
        vm = utils.ValuesMap.from_timeseries(arr, nnodes=nnodes, nticks=nticks)
        self.assertEqual(vm.shape, (nticks, nnodes))
        # The pattern should repeat 5 times
        for i in range(nticks):
            expected_value = arr[i % 10]
            self.assertTrue(np.all(vm.values[i] == expected_value))

    def test_from_timeseries_replicate_partial(self):
        """
        Validate ValuesMap.from_timeseries() handles partial pattern repetition correctly.

        SCENARIO SETUP:
        - Time series: [0, 1, 2, ..., 9] (10-tick pattern)
        - Requested nticks: 25 (2 full cycles + 5 extra ticks)
        - Nodes: 3

        WHAT IS TESTED:
        - Output shape is (25, 3) as requested
        - Pattern repeats 2 full times (ticks 0-19) then 5 additional values (ticks 20-24)
        - Modulo arithmetic correctly handles incomplete final cycle
        - No discontinuity or padding at the truncation point

        FAILURE MEANING:
        If this test fails, from_timeseries() isn't correctly handling partial repetitions.
        This would break scenarios where simulation length isn't an exact multiple of the
        seasonal pattern length (e.g., 912-day (2.5 year) simulation with 365-day yearly pattern).
        Incorrect handling could cause abrupt jumps or wrong values at cycle boundaries.
        """
        arr = np.arange(10, dtype=np.float32)
        nticks = 25
        nnodes = 3
        vm = utils.ValuesMap.from_timeseries(arr, nnodes=nnodes, nticks=nticks)
        self.assertEqual(vm.shape, (nticks, nnodes))
        # The pattern should repeat 2 full times and then 5 more values
        for i in range(nticks):
            expected_value = arr[i % 10]
            self.assertTrue(np.all(vm.values[i] == expected_value))


class TestTimingStats(unittest.TestCase):
    """
    Test suite for _TimingStats utility that profiles code execution timing.

    _TimingStats provides hierarchical performance profiling via context managers.
    Useful for identifying performance bottlenecks in model components and simulations.

    WHAT IS TESTED:
    - Context manager syntax (with stats.start("label"))
    - Freezing stats to finalize timing measurements
    - String output format for human-readable timing reports
    - Dictionary output format for programmatic access to timing data
    """

    def test_timing_stats(self):
        """
        Validate _TimingStats basic functionality for timing code blocks.

        SCENARIO SETUP:
        - Create TimingStats object
        - Time a simple loop block labeled "test"
        - Freeze stats to finalize measurements
        - Extract results as both string and dictionary

        WHAT IS TESTED:
        - Context manager correctly times code execution
        - freeze() successfully finalizes statistics
        - to_string() produces output containing the label "test"
        - to_dict() produces structured dict with "children" containing labeled timing node

        FAILURE MEANING:
        If this test fails, the timing infrastructure is broken. String output failure suggests
        formatting issues. Dictionary output failure (missing label in children) indicates the
        timing tree structure isn't being built correctly. This would prevent performance
        profiling and optimization of model components.
        """
        stats = utils._TimingStats()
        with stats.start("test"):
            time_sum = 0
            for _ in range(10):
                time_sum += 1
        stats.freeze()
        s = stats.to_string()
        self.assertIn("test", s)
        d = stats.to_dict()
        self.assertEqual(d["children"][0]["label"], "test")


class TestValidateDecorator(unittest.TestCase):
    """
    Test suite for @validate decorator that wraps methods with pre/post validation hooks.

    The @validate decorator enables runtime validation by calling pre-condition and
    post-condition functions before and after a method executes. Used to ensure
    model consistency and catch state errors during simulation.

    WHAT IS TESTED:
    - Decorator correctly calls pre-hook before method
    - Decorator correctly calls post-hook after method
    - Execution order: pre → method → post
    - Decorator only activates when model.validating = True
    """

    def test_validate_decorator(self):
        """
        Validate @validate decorator invokes pre/post hooks in correct order.

        SCENARIO SETUP:
        - Mock class with validating=True to enable validation
        - Method decorated with @validate(pre, post)
        - Tracking list captures function call order

        WHAT IS TESTED:
        - pre() hook called before method execution
        - post() hook called after method execution
        - Execution sequence is exactly: ["pre5", "run5", "post5"]
        - Argument (tick=5) correctly passed to all three functions

        FAILURE MEANING:
        If this test fails, the @validate decorator isn't wrapping methods correctly. Wrong
        execution order suggests pre/post hooks aren't being called at the right time. Missing
        calls suggest the decorator isn't invoking hooks at all. This would break validation
        infrastructure, allowing invalid model states to persist undetected during simulations,
        potentially causing subtle bugs or incorrect results.
        """
        calls = []

        class Dummy:
            def __init__(self):
                self.model = self
                self.validating = True

            def pre(self, tick):
                calls.append(f"pre{tick}")

            def post(self, tick):
                calls.append(f"post{tick}")

            @utils.validate(pre, post)
            def run(self, tick):
                calls.append(f"run{tick}")

        d = Dummy()
        d.run(5)
        self.assertEqual(calls, ["pre5", "run5", "post5"])

    class TestValidateDecoratorModelVsSelf(unittest.TestCase):
        """
        Additional tests for @validate decorator focusing on model vs self validation flags.

        Specifically, test that validation occurs when the object's model has validating=True,
        even if the object itself does not have validating=True.
        """

        def test_validate_decorator_model_validating(self):
            """
            Validate @validate decorator triggers validation when model.validating=True,
            even if self.validating=False.

            SCENARIO SETUP:
            - Dummy object with self.validating=False, but self.model.validating=True
            - Method decorated with @validate(pre, post)
            - Tracking list captures function call order

            WHAT IS TESTED:
            - pre() hook called before method execution
            - post() hook called after method execution
            - Execution sequence is exactly: ["pre7", "run7", "post7"]
            - Argument (tick=7) correctly passed to all three functions

            FAILURE MEANING:
            If this test fails, the @validate decorator is not correctly checking model.validating.
            This would break validation in scenarios where validation is controlled at the model
            level rather than per-component.
            """
            calls = []

            class Dummy:
                class Model:
                    def __init__(self):
                        self.validating = True

                def __init__(self):
                    self.model = Dummy.Model()
                    self.validating = False  # self is not validating, but model will be set to True

                def pre(self, tick):
                    calls.append(f"pre{tick}")

                def post(self, tick):
                    calls.append(f"post{tick}")

                @utils.validate(pre, post)
                def run(self, tick):
                    calls.append(f"run{tick}")

            d = Dummy()
            d.model.validating = True  # Set model.validating to True
            d.run(7)
            self.assertEqual(calls, ["pre7", "run7", "post7"])


class TestGetCentroids(unittest.TestCase):
    """
    Test suite for get_centroids() function that computes polygon centroids.

    get_centroids() extracts centroid points from GeoDataFrame polygon geometries,
    used for calculating distances between spatial nodes in gravity models and
    for spatial visualization.

    WHAT IS TESTED:
    - Centroid extraction from single polygon
    - Centroid extraction from multiple polygons
    - CRS (coordinate reference system) preservation
    """

    def test_get_centroids_multiple(self):
        """
        Validate get_centroids() correctly processes multiple polygons.

        SCENARIO SETUP:
        - Two identical unit square polygons: [(0,0), (1,0), (1,1), (0,1)]
        - GeoDataFrame with CRS = EPSG:4326 (WGS84 lat/lon)

        WHAT IS TESTED:
        - Output contains exactly 2 centroid points (one per polygon)
        - CRS preserved correctly (output CRS = 4326)
        - Implicit: centroids calculated at (0.5, 0.5) for unit square

        FAILURE MEANING:
        If this test fails, get_centroids() isn't correctly extracting centroids from
        multi-polygon GeoDataFrames. Count mismatch suggests polygon iteration is broken.
        CRS mismatch suggests coordinate system isn't being preserved, which would break
        distance calculations in gravity models. Spatial transmission would use wrong
        distances, producing incorrect epidemic spread patterns.
        """
        poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        gdf = gpd.GeoDataFrame({"geometry": [poly, poly]}, crs="EPSG:4326")
        centroids = utils.get_centroids(gdf)
        self.assertEqual(len(centroids), 2)
        self.assertEqual(centroids.crs.to_epsg(), 4326)

    def test_get_centroids_single(self):
        """
        Validate get_centroids() correctly processes a single polygon.

        SCENARIO SETUP:
        - One unit square polygon: [(0,0), (1,0), (1,1), (0,1)]
        - GeoDataFrame with CRS = EPSG:4326 (WGS84 lat/lon)

        WHAT IS TESTED:
        - Output contains exactly 1 centroid point
        - CRS preserved correctly (output CRS = 4326)
        - Edge case: single-polygon input doesn't cause errors

        FAILURE MEANING:
        If this test fails, get_centroids() fails on single-polygon input. This is an
        important edge case since single-node models are common for testing and simple
        scenarios. Failure suggests the function assumes multiple polygons or has issues
        with array reshaping/iteration. CRS failure has same implications as multi-polygon test.
        """
        poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        gdf = gpd.GeoDataFrame({"geometry": [poly]}, crs="EPSG:4326")
        centroids = utils.get_centroids(gdf)
        self.assertEqual(len(centroids), 1)
        self.assertEqual(centroids.crs.to_epsg(), 4326)


class TestGetDefaultParameters(unittest.TestCase):
    """
    Test suite for get_default_parameters() function that provides baseline parameter sets.

    get_default_parameters() returns a dictionary of standard parameter values for
    epidemic simulations, providing sensible defaults to reduce boilerplate in model setup.

    WHAT IS TESTED:
    - Presence of required parameter keys
    - Correct default values for critical parameters
    - Parameter types and structure
    """

    def test_get_default_parameters(self):
        """
        Validate get_default_parameters() returns correct default parameter dictionary.

        WHAT IS TESTED:
        - Dictionary contains required keys: nticks, beta, verbose
        - Default values match expected baseline configuration:
          - nticks = 730 (2 years for observing seasonal patterns)
          - beta = 0.15 (moderate transmission rate)
          - verbose = False (quiet output by default)

        FAILURE MEANING:
        If this test fails, default parameters have changed or are missing. This breaks
        backward compatibility and could cause existing code to fail. Missing keys suggest
        incomplete parameter set. Wrong values suggest defaults were changed without updating
        tests. Since many examples and tutorials rely on these defaults, changes need to be
        carefully managed. Failure indicates need to update documentation and example code
        if the defaults have been intentionally changed.
        """
        params = utils.get_default_parameters()
        self.assertIn("nticks", params)
        self.assertIn("beta", params)
        self.assertEqual(params["nticks"], 730)
        self.assertEqual(params["beta"], 0.15)
        self.assertEqual(params["verbose"], False)


class TestSeedingFunctions(unittest.TestCase):
    def setUp(self):
        self.model = DummyModel(10)

    def test_seed_infections_randomly(self):
        """
        Test that seed_infections_randomly() correctly infects a specified number of agents
        across the entire people.

        Test design:
        - Start with a people of 10 fully susceptible agents.
        - Request 5 random infections.
        - Validate that exactly 5 agents transition to INFECTIOUS.
        - Confirm that their itimers are initialized to the model's inf_mean.

        Pass = The correct number of agents are infected and initialized properly.
        Fail = Too few/many infected, incorrect states, or incorrect timer values.
        """
        nodeids = utils.seed_infections_randomly(self.model, ninfections=5)
        self.assertEqual(len(nodeids), 5)
        self.assertEqual(np.sum(self.model.people.state == State.INFECTIOUS.value), 5)
        self.assertTrue(np.all(self.model.people.itimer[nodeids] == self.model.params.inf_mean))

    def test_seed_infections_in_patch(self):
        """
        Test that seed_infections_in_patch() correctly infects a specified number of agents
        within a single node (patch), using the state array rather than susceptibility.

        Test design:
        - Use a people of 10 agents assigned to 5 nodes (2 agents per node).
        - Request 2 infections in node 2 (ipatch=2).
        - Validate that exactly 2 individuals with nodeid==2 transition to
          INFECTIOUS state (state == State.INFECTIOUS.value).
        - Confirm that their itimers are initialized to inf_mean.
        - Optionally verify that no agents in other nodes become infectious.

        Pass =
        - Exactly 2 agents in node 2 are marked INFECTIOUS and receive correct timers.
        - No agents in other nodes are marked INFECTIOUS.

        Fail =
        - Infections spill to other nodes,
        - Wrong number of agents in node 2 become infectious,
        - Or their timers are not initialized to the expected value.
        """
        # Arrange: 5 nodes, 2 agents per node
        self.model.people.nodeid = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])

        # Act: seed infections only in patch 2
        utils.seed_infections_in_patch(self.model, ipatch=2, ninfections=2)

        # Agents in node 2 that are INFECTIOUS
        infected_in_patch = (self.model.people.nodeid == 2) & (self.model.people.state == State.INFECTIOUS.value)

        # Agents in other nodes that are INFECTIOUS (should be none)
        infected_elsewhere = (self.model.people.nodeid != 2) & (self.model.people.state == State.INFECTIOUS.value)

        # Assert: exactly two infectious in node 2
        self.assertEqual(np.sum(infected_in_patch), 2)

        # Assert: no infectious agents outside node 2
        self.assertEqual(np.sum(infected_elsewhere), 0)

        # Assert: their timers are set to inf_mean
        self.assertTrue(np.all(self.model.people.itimer[infected_in_patch] == self.model.params.inf_mean))


if __name__ == "__main__":
    unittest.main()
