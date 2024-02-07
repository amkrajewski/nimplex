# Copyrigth (C) 2023 Adam M. Krajewski
# License: MIT

# Pragmas with compiler and linker options
{.passC: "-flto -ffast-math".} 
{.passL: "-flto".} 

# Standard Nim library imports
from std/math import binom, ln
import std/sugar
import std/times
import std/strutils
from std/algorithm import reverse
from std/sequtils import foldl, mapIt

# Arraymancer library for tensor operations
import arraymancer/Tensor
import arraymancer/io

# OS module can cause issues for Python bindings, so it's not imported in the library mode, when it is not needed.
when appType != "lib":
    import std/os

# Nimpy module for Python bindings when running in the library mode and not generating documentation (to avoid duplicate API entries)
when appType == "lib" and not defined(nimdoc):
    import nimpy

when defined(nimdoc):
    # All of (comprehensive) introduction to the documentation lives in this included Nim file, while API is generated from docstrings in the code. It was moved there for cleaner code.
    include docs/docs
    # The plotting utils are not part of the core library, but are imported during documentation generation to index them as part of the library.
    import utils/plotting

# GRID
proc simplex_grid*(dim: int, 
                   ndiv: int): Tensor[int] =
    ## .. image:: ../assets/small_FI.png               
    ## Generates a full (including the simplex boundary) simplex grid in a `dim`-component space with `ndiv` divisions per dimension (i.e., quantized to `1/ndiv`).
    ## The result is a deterministically allocated Arraymancer `Tensor[int]` of shape `(N_S(dim, ndiv), dim)` containing all possible compositions in the simplex space
    ## expressed as integer numbers of quantum units. The grid is generated procedurally using a modified version of NEXCOM algorithm (see manuscript for details).
    # L is the total number of unique points in the simplex grid, which we know a priori
    let L: int = binom(ndiv+dim-1, dim-1)
    result = newTensor[int]([L, dim])
    # x is the current composition
    var x = zeros[int](dim)
    x[dim-1] = ndiv
    for j in 0..<dim:
        result[0, j] = x[j]
    var h = dim
    for i in 1..<L:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..<dim:
            result[i, j] = x[j]
        if val != 1:
            h = dim
    return result

proc simplex_grid_fractional*(dim: int,
                              ndiv: int): Tensor[float] =
    ## .. image:: ../assets/small_FF.png
    ## Conceptually similar to `simplex_grid`_ but results in
    ## Arraymancer `Tensor[float]` of shape `(N_S(dim, ndiv), dim)` containing all possible compositions in the simplex space
    ## expressed as fractions summing to 1. The grid is generated procedurally using a modified version of NEXCOM algorithm (see manuscript for details) and then normalized
    ## to 1 in a quick single-pass post-processing step.
    result = simplex_grid(dim, ndiv).asType(float)
    result = result.map(x => x / float(ndiv))
    return result


proc simplex_internal_grid*(dim: int, 
                            ndiv: int): Tensor[int] =
    ## Same as `simplex_grid`_ but generates only points inside the simplex, i.e., all components are present.
    # L is the total number of unique points inside the simplex grid, which we know a priori
    let L: int = binom(ndiv-1, dim-1)
    result = newTensor[int]([L, dim])
    var x = ones[int](dim)
    x[dim-1] = ndiv+1-dim
    for j in 0..dim-1:
        result[0, j] = x[j]
    var h = dim
    for i in 1..<L:
        h -= 1
        let val = x[h]
        x[h] = 1
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..dim-1:
            result[i, j] = x[j]
        if val != 2:
            h = dim
    return result

proc simplex_internal_grid_fractional*(dim: int,
                                       ndiv: int): Tensor[float] =
    ## Conceptually similar to `simplex_internal_grid`_ but results in 
    ## Arraymancer `Tensor[float]` of shape `(N_I(dim, ndiv), dim)` containing all possible compositions in the simplex space
    ## expressed as fractions summing to 1 (same as ``simplex_grid_fractional``). 
    result = simplex_internal_grid(dim, ndiv).asType(float)
    result = result.map(x => x / float(ndiv))
    return result

# RANDOM SAMPLING

proc simplex_sampling_mc*(dim: int,
                         samples: int): Tensor[float] =
    ## Randomly samples `samples` number of points from a simplex in a `dim`-component space. The result is a deterministically allocated Arraymancer `Tensor[float]` of shape `(samples, dim)` containing all sampled points
    ## expressed as fractions summing to 1. As eluded to in [Capabilities](#capabilities), this is not as straightforward as in Euclidean spaces.
    ## The sampling is done by generating `samples` number of points from a special case of Dirichlet distribution (alpha=1) by sampling from a uniform distribution in (dim)-Cartesian space, then
    ## taking the negative logarithm of each coordinate, and finally normalizing each point to sum to 1. This approach is described in the manuscript and contrasted with another (lower performance) approach involving sorting.
    let neglograndom = randomTensor[float](
        [samples, dim], 
        1.0
        ).map(x => -ln(x))
    let sums = neglograndom.sum(axis=1)
    result = neglograndom /. sums

# GRAPH

proc simplex_graph_3C*(
    ndiv: int): (Tensor[int], seq[seq[int]]) =
    ## A special case of `simplex_graph`_ for 3-component spaces, which is optimized for performance by avoiding the use of `binom` function and filling the neighbor list sequentially only once per node (in both directions).
    let L: int = binom(ndiv+2, 2)
    var nodes = newTensor[int]([L, 3])
    var neighbors = newSeq[seq[int]](L)
    var x = zeros[int](3)

    func neighborsLink(i:int, x:Tensor, ndiv:int): seq[int] =
        const jump0 = 1
        let jump1 = 1+ndiv-x[0]

        if x[0] != 0:
            result.add(i-jump1)
            result.add(i-jump1-jump0)
        if x[1] != 0:
            result.add(i-jump0)
            result.add(i+jump1-jump0)
        if x[2] != 0:
            result.add(i+jump0)
            result.add(i+jump1)
        return result

    x[2] = ndiv
    for j in 0..2:
        nodes[0, j] = x[j]
    var h = 3

    neighbors[0] = neighborsLink(0, x, ndiv)

    for i in 1..<L:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[2] = val - 1
        x[h-1] += 1
        for j in 0..2:
            nodes[i, j] = x[j]
        neighbors[i] = neighborsLink(i, x, ndiv)
        if val != 1:
            h = 3
    return (nodes, neighbors)

proc simplex_graph_3C_fractional*(ndiv: int): (Tensor[float], seq[seq[int]]) =
    ## A special case of `simplex_graph_fractional`_ for 3-component spaces utilizing optimized `simplex_graph_3C`_.
    let graph = simplex_graph_3C(ndiv)
    var nodes = graph[0].asType(float)
    nodes = nodes.map(x => x / float(ndiv))
    return (nodes, graph[1])

proc simplex_graph*(
    dim: int, 
    ndiv: int): (Tensor[int], seq[seq[int]]) =
    ## .. image:: ../assets/small_GI.png
    ## Generates a simplex graph in a `dim`-component space based on (1) grid of nodes following `ndiv` divisions per dimension (i.e., quantized to `1/ndiv`),
    ## similar to `simplex_grid`_, and (2) a list of neighbor lists corresponding to edges. The result is a tuple of 
    ## (1) a deterministically allocated Arraymancer `Tensor[int]` of shape `(N_S(dim, ndiv), dim)` containing all possible compositions in the simplex space
    ## just like in `simplex_grid`_, and (2) a `seq[seq[int]]` containing a list of neighbors for each node. The current implementation utilizes GC-allocated `seq` for neighbors
    ## to reduce memory footprint in cases where `ndiv` is close to `dim` and not all nodes have the complete set of `dim(dim-1)` neighbors. This is a tradeoff between memory and performance, which can be
    ## adjusted by switching to a (N_S(dim, ndiv), dim(dim-1)) `Tensor[int]` with `-1` padding for missing neighbors, just like done in `outFunction_graph`_ for NumPy output generation.
    let L: int = binom(ndiv+dim-1, dim-1)
    var nodes = newTensor[int]([L, dim])
    var neighbors = newSeq[seq[int]](L)
    var x = zeros[int](dim)

    proc neighborsLink(i:int, x:Tensor, ndiv:int, dim:int, 
                       neighbors: var seq[seq[int]]): void =
        var jumps = newSeq[int](dim-1)
        jumps[0] = 1  #binom(x,0)=1
        for j in 1..<(dim-1):
            jumps[j] = binom(j+ndiv-sum(x[0..(dim-2-j)]), j)
        var temp: int
        for order in 0..(dim-2): 
            temp = 0
            if x[order] != 0:
                for dir in 0..(dim-2-order): 
                    temp += jumps[dim-2-order-dir]
                    neighbors[i].add(i - temp)
                    neighbors[i - temp].add(i)

    x[dim-1] = ndiv
    for j in 0..<dim:
        nodes[0, j] = x[j]
    var h = dim

    neighborsLink(0, x, ndiv, dim, neighbors)

    for i in 1..<L:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..<dim:
            nodes[i, j] = x[j]
        neighborsLink(i, x, ndiv, dim, neighbors)
        if val != 1:
            h = dim
    return (nodes, neighbors)

proc simplex_graph_fractional*(dim: int, ndiv: int): (Tensor[float], seq[seq[int]]) =
    ## Conceptually similar to `simplex_graph`_ but the first part of the result tuple (graph nodes) is an Arraymancer `Tensor[float]` of shape `(N_S(dim, ndiv), dim)` containing all possible compositions in the simplex space
    ## expressed as fractions summing to 1 (same as `simplex_grid_fractional`).
    let graph = simplex_graph(dim, ndiv)
    var nodes = graph[0].asType(float)
    nodes = nodes.map(x => x / float(ndiv))
    return (nodes, graph[1])

# CORE UTILS

proc attainable2elemental*(simplexPoints: Tensor[float],
                           components: seq[seq[float]]): Tensor[float] =
    ## Accepts a `simplexPoints` Arraymancer `Tensor[float]` of shape corresponding to a simplex grid (e.g., from `simplex_grid_fractional`_) or random samples (e.g., from `simplex_sampling_mc`_) and a `components` list of lists of floats, which represents a list of
    ## compositions in the **elemental** space serving as base components of the **attainable** space given in `simplexPoints`. The `components` can be a row-consistnet mixed list list of integer and fractional compositions, to allow for both types of inputs. 
    ## It then projects each point from the attainable space to the elemental space using matrix multiplication.
    runnableExamples:
        const components = @[
            @[0.94, 0.05, 0.01], # Fe95 C5 Mo1
            @[3.0, 1.0, 0.0],    # Fe3C
            @[0.2, 0.0, 0.8]     # Fe20 Mo80
        ]
        let grid = simplex_grid_fractional(3, 4)
        let elementalGrid = grid.attainable2elemental(components)
        echo elementalGrid
    # Tensor of components which can be "integer" ([2,2,1]) or "fractional" ([0.4,0.4,0.2]) compositions
    var cmpTensor: Tensor[float] = components.toTensor()
    # Normalize components to sum to 1
    cmpTensor = cmpTensor /. cmpTensor.sum(axis=1)
    # Matrix multiplication to get the final grid
    result = simplexPoints * cmpTensor

func pure_component_indexes*(dim: int, ndiv: int): seq[int] =
    ## This helper function returns a `seq[int]` of indexes of pure components in a simplex grid of `dim` dimensions and `ndiv` divisions per dimension (e.g., from `simplex_grid`_).
    for d in 1..dim:
        result.add(binom(ndiv+d-1, ndiv)-1)
    # Reverse the order as the last pure component is the first in the grid
    result.reverse()

func pure_component_indexes_internal*(dim: int, ndiv: int): seq[int] =
    ## This helper function returns a `seq[int]` of indexes of pure components in an **internal** simplex grid of `dim` dimensions and `ndiv` divisions per dimension (e.g., from `simplex_internal_grid`_).
    for d in 1..dim:
        result.add(binom(ndiv-1, dim-1)-1)
    # Reverse the order as the last pure component is the first in the grid
    result.reverse()


# PYTHON BINDINGS
when appType == "lib" and not defined(nimdoc):
    # Direct translation of the Nim API to Python using Nimpy
    
    proc simplex_grid_py*(dim: int, ndiv: int): seq[seq[int]] {.exportpy.} = 
        simplex_grid(dim, ndiv).toSeq2D()

    proc simplex_grid_fractional_py*(dim: int, ndiv: int): seq[seq[float]] {.exportpy.} = 
        simplex_grid_fractional(dim, ndiv).toSeq2D()

    proc simplex_internal_grid_py*(dim: int, ndiv: int): seq[seq[int]] {.exportpy.} = 
        simplex_internal_grid(dim, ndiv).toSeq2D()

    proc simplex_internal_grid_fractional_py*(dim: int, ndiv: int): seq[seq[float]] {.exportpy.} = 
        simplex_internal_grid_fractional(dim, ndiv).toSeq2D()

    proc simplex_sampling_mc_py*(dim: int, samples: int): seq[seq[float]] {.exportpy.} = 
        simplex_sampling_mc(dim, samples).toSeq2D() 

    proc simplex_graph_3C_py*(ndiv: int): (seq[seq[int]], seq[seq[int]]) {.exportpy.} = 
        let graph = simplex_graph_3C(ndiv)
        return (graph[0].toSeq2D(), graph[1])

    proc simplex_graph_3C_fractional_py*(ndiv: int): (seq[seq[float]], seq[seq[int]]) {.exportpy.} =
        let graph = simplex_graph_3C_fractional(ndiv)
        return (graph[0].toSeq2D(), graph[1])

    proc simplex_graph_py*(dim: int, ndiv: int): (seq[seq[int]], seq[seq[int]]) {.exportpy.} = 
        let graph = simplex_graph(dim, ndiv)
        return (graph[0].toSeq2D(), graph[1])

    proc simplex_graph_fractional_py*(dim: int, ndiv: int): (seq[seq[float]], seq[seq[int]]) {.exportpy.} =
        let graph = simplex_graph_fractional(dim, ndiv)
        return (graph[0].toSeq2D(), graph[1])

    # Extra functions for Python bindings only
    proc embeddedpair_simplex_grid_fractional_py*(components: seq[seq[float]], ndiv: int): (seq[seq[float]], seq[seq[float]]) {.exportpy.} =
        let temp = simplex_grid_fractional(components.len, ndiv)
        result[0] = temp.toSeq2D()
        result[1] = temp.attainable2elemental(components).toSeq2D()

    proc embeddedpair_simplex_internal_grid_fractional_py*(components: seq[seq[float]], ndiv: int): (seq[seq[float]], seq[seq[float]]) {.exportpy.} =
        let temp = simplex_internal_grid_fractional(components.len, ndiv)
        result[0] = temp.toSeq2D()
        result[1] = temp.attainable2elemental(components).toSeq2D()

    proc embeddedpair_simplex_sampling_mc_py*(components: seq[seq[float]], samples: int): (seq[seq[float]], seq[seq[float]]) {.exportpy.} =
        let temp = simplex_sampling_mc(components.len, samples)
        result[0] = temp.toSeq2D()
        result[1] = temp.attainable2elemental(components).toSeq2D()

    proc embeddedpair_simplex_graph_3C_fractional_py*(components: seq[seq[float]], ndiv: int): (seq[seq[float]], seq[seq[float]], seq[seq[int]]) {.exportpy.} =
        let temp = simplex_graph_3C_fractional(ndiv)
        result[0] = temp[0].toSeq2D()
        result[1] = temp[0].attainable2elemental(components).toSeq2D()
        result[2] = temp[1]

    proc embeddedpair_simplex_graph_fractional_py*(components: seq[seq[float]], ndiv: int): (seq[seq[float]], seq[seq[float]], seq[seq[int]]) {.exportpy.} =
        let temp = simplex_graph_fractional(components.len, ndiv)
        result[0] = temp[0].toSeq2D()
        result[1] = temp[0].attainable2elemental(components).toSeq2D()
        result[2] = temp[1]

    proc pure_component_indexes_py*(dim: int, ndiv: int): seq[int] {.exportpy.} =
        pure_component_indexes(dim, ndiv)

    proc pure_component_indexes_internal_py*(dim: int, ndiv: int): seq[int] {.exportpy.} =
        pure_component_indexes_internal(dim, ndiv)


# UTILS

template benchmark(benchmarkName: string, code: untyped) =
    ## A simple benchmarking template which takes a name and a code block to run. It prints the benchmark name and the duration of the code block execution in natural time units
    ## quantized to microseconds (e.g., 119 milliseconds and 296 microseconds). All tasks together should take 500-1000 milliseconds on a modern CPU.
    block:
        let t0 = cpuTime()
        code
        let t1 = cpuTime()
        echo benchmarkName & "\n" & $initDuration(microseconds = ((t1 - t0)*1e6).int) & "\n"

proc echoHelp*() = 
    ## Prints the help message for the CLI, which is a concise version of one given in nimplex's documentation.
    echo """

To run nimplex please either (1) provide no arguments and follow the prompts or 
(2) use "-c" or "--config" to provide the configuration per instructions below:

- Provide the 3-letter configuration for task type:
    1. Grid type or uniform random sampling:
        - F: Full grid (including the simplex boundary)
        - I: Internal grid (only points inside the simplex)
        - R: Random/Monte Carlo uniform sampling over simplex.
        - G: Graph (list of grid nodes and list of their neighbors)
    2. Fractional or Integer positions:
        - F: Fractional grid/graph (points are normalized to fractions of 1)
        - I: Integer grid/graph (points are integers)
    3. Print full result, its shape, or persist in a file:
        - P: Print (presents full result as a table)
        - S: Shape (only the shape / size information)
        - N: Persist to NumPy array file ("nimplex_<configFlags>.npy" or 
             optionally a custom path as an additonal argument)

- Followed by integers of (1) simplex dimension and (2) number of divisions or
  samples depending on the task type. Optionally, custom output file path for 
  NumPy array can be provided as the last argument. E.g.:
    -c FFS [simplex dimension] [number of divisions]
    -c RFP [simplex dimension] [number of samples]
    -c FIN [simplex dimension] [number of divisions] [path/to/outfile.npy]

You can also utilize the following auxiliary flags:
--help       | -h   --> Show help.
--benchmark  | -b   --> Run benchmark for all tasks (9-dimensional space
                        with 12 divisions per dimension / 1M random samples).
"""

func configValidation(config: string) = 
    ## Validates the 3-letter configuration string provided by the user.
    assert config.len == 3, "\n--> Invalid configuration lenght. Must be 3 letters."
    assert config[0] in @['F', 'I', 'R', 'G'], "\n--> Invalid configuration (in the 1st letter). Must be F, I or R for Full grid, Internal grid, Monte Carlo sampling, or Graph respectively"
    assert config[1] in @['F', 'I'], "\n--> Invalid configuration (in the 2nd letter). Must be F, or I for Fractional positions, or Integer positions respectively"
    if config[0] == 'R':
        assert config[1] == 'F', "\n--> Integer positions not implemented for Random sampling. Must be F for Fractional positions."
    assert config[2] in @['P', 'S', 'N'], "\n--> Invalid configuration (in the 3rd letter). Must be P, S or N for Print full result, Shape, or persist Numpy output respectively"

func nDivValidation(config: string, nDiv: int, dim: int) = 
    ## Validates the number of divisions per each simplex dimension provided by the user for all tasks except Random sampling.
    if config[0] == 'I':
        assert ndiv >= dim, "\n--> Invalid number of divisions. Must be greater or equal to the simplex dimension to produce a non-empty internal grid."
    else:
        assert ndiv > 0, "\n--> Invalid number of divisions. Must be a positive integer."

proc outFunction(config: string, npyName: string, outputData: Tensor) =
    ## Handles the output of a simplex grid or random sampling in `outputData` `Tensor` when run from the CLI, based on the 3rd letter of the configuration string `config`,
    ## and the destiantion filename `npyName`, which should include the extension.
    case config[2]:
        of 'P': 
            echo "Full Output: ", outputData
        of 'N': 
            outputData.write_npy(npyName)
            echo "Shape: ", outputData.shape
        of 'S': 
            echo "Shape: ", outputData.shape
        else: 
            echo "Invalid Config"
    

proc outFunction_graph(config: string, dim: int, ndiv: int, npyName: string, outputData: (Tensor, seq)) =
    ## Handles the output of graph tasks when run from the CLI, based on the 3rd letter of the configuration string `config`, the number of dimensions `dim`, the number of divisions per dimension `ndiv`, which
    ## together are used to determine the graph edges Tensor shape (using math) faster than querying the input edges sequence. The destiantion filename `npyName`, which should include the extension, also needs to be provided.
    case config[2]:
        of 'P': 
            echo "Nodes:"
            echo outputData[0]
            echo "Neighbors:"
            echo outputData[1]
        of 'N': 
            outputData[0].write_npy(npyName.replace(".npy", "_nodes.npy"))
            let 
                maxNeighbors = dim*(dim-1)
                L = binom(ndiv+dim-1, dim-1)
            var neighborsTensor = newTensor[int]([L, maxNeighbors])
            for i in 0..<L:
                for j in 0..<maxNeighbors:
                    if j < outputData[1][i].len:
                        neighborsTensor[i, j] = outputData[1][i][j]
                    else:
                        neighborsTensor[i, j] = -1
            neighborsTensor.write_npy(npyName.replace(".npy", "_neighbors.npy"))
            echo "Nodes Shape: ", outputData[0].shape
            echo "Edges Count: ", outputData[1].mapIt(it.len).foldl(a + b, 0)
        of 'S': 
            echo "Nodes Shape: ", outputData[0].shape
            echo "Edges Count: ", outputData[1].mapIt(it.len).foldl(a + b, 0)
        else: 
            echo "Invalid Congig"

proc taskRouter(config: string, dim: int, ndiv: int, npyName: string) =
    ## Routes the task to the appropriate calculation and output function based on the first 2 letters of the configuration string.
    case config[0..1]:
        of "FF": outFunction(
            config, npyName, simplex_grid_fractional(dim, ndiv))
        of "FI": outFunction(
            config, npyName, simplex_grid(dim, ndiv))
        of "IF": outFunction(
            config, npyName, simplex_internal_grid_fractional(dim, ndiv))
        of "II": outFunction(
            config, npyName, simplex_internal_grid(dim, ndiv))
        of "RF": outFunction(
            config, npyName, simplex_sampling_mc(dim, samples=ndiv))
        of "GI": outFunction_graph(
            config, dim, ndiv, npyName, simplex_graph(dim, ndiv))
        of "GF": outFunction_graph(
            config, dim, ndiv, npyName, simplex_graph_fractional(dim, ndiv))
        else:
            echo "\n--> Invalid configuration in the first 2 config letters."
            quit(1)

when appType != "lib":
    when isMainModule:
        let args = commandLineParams() ## \
        ## Command line arguments parsed when the module is run as a script, rather than as a library, allowing efortless CLI usage without any Python or Nim knowledge.
        ## When empty, interactive mode is triggered and user is navigated through the configuration process. Otherwise, the first argument is expected to be `-c` or `--config` 
        ## followed by the configuration flags and parameters as described in the help message below. See `echoHelp()` for more details.
        
        # Interactive
        if args.len == 0:
            echo "Configuration (Full/Internal/Random/Graph)(Fractional/Integer)(Print/Shape/Numpysave) - e.g. FFS/RFP/FIN:"
            let config = readLine(stdin)
            configValidation(config)

            echo "Simplex Dimensions / N of Components:"
            let dim = readLine(stdin).parseInt() 

            var nDiv: int
            if config[0]=='R':
                echo "Number of Samples:"
                nDiv = readLine(stdin).parseInt()
                assert nDiv > 0, "\n--> Invalid number of samples. Must be a positive integer"
            else:
                echo "N Divisions per Dimension:"
                ndiv = readLine(stdin).parseInt() 
                nDivValidation(config, ndiv, dim)

            var npyName: string = "nimplex_" & config[0..1] & "_" & $dim & "_" & $ndiv & ".npy"
            if config[2] == 'N':
                echo "NumPy Array Output Filename (skip for default: " & npyName & "):"
                let tempIn = readLine(stdin)
                if tempIn.len > 0:
                    npyName = tempIn
                echo "Persisting to NumPy array file: ", npyName

            taskRouter(config, dim, ndiv, npyName)

        # Configured
        elif args[0] == "-c" or args[0] == "--config":
            let config = args[1]
            echo "Running with configuration: ", args[1..<args.len]
            configValidation(config)

            let dim = args[2].parseInt()
            assert dim > 0, "Invalid dimension"

            var nDiv: int
            if config[0]=='R':
                nDiv = args[3].parseInt()
                assert nDiv > 0, "\n--> Invalid sample number. Must be a positive integer"
            else:
                ndiv = args[3].parseInt()
                nDivValidation(config, ndiv, dim)

            var npyName: string = "nimplex_" & config[0..1] & "_" & $dim & "_" & $ndiv & ".npy"
            if config[2] == 'N':
                if args.len == 5:
                    npyName = args[4]
                echo "Persisting to NumPy array file: ", npyName

            taskRouter(config, dim, ndiv, npyName)

        elif args[0] in @["-h", "--help"]:
            echoHelp()
            quit(0)

        elif args[0] in @["-b", "--benchmark"]:
            # A few benchmarks for the library to compare across different systems, implementations, and languages.
            benchmark "Simplex Grid Full (dim=9, ndiv=12):":
                discard simplex_grid(9, 12)
            benchmark "Simplex Grid Full Fractional (dim=9, ndiv=12):":
                discard simplex_grid_fractional(9, 12)
            benchmark "Simplex Grid Internal (dim=9, ndiv=12):":
                discard simplex_internal_grid(9, 12)
            benchmark "Simplex Random Sampling (dim=9, samples=1M):":
                discard simplex_sampling_mc(9, 1_000_000)
            benchmark "Simplex Graph (dim=9, ndiv=12):":
                discard simplex_graph(9, 12)
            benchmark "Simplex Graph ND (dim=3, ndiv=1000):":
                discard simplex_graph(3, 1000)
            benchmark "Simplex Graph 3D (dim=3, ndiv=1000):":
                discard simplex_graph_3C(1000)

        # Fallback
        else:
            echoHelp()
            quit(1)