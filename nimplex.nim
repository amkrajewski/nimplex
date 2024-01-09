# Copyrigth (C) 2023 Adam M. Krajewski

from std/math import binom, ln
import std/sugar
import std/os
import std/times
import std/strutils

import arraymancer/Tensor
import arraymancer/io

import nimpy

proc simplex_grid*(dim: int, 
                   ndiv: int): Tensor[int] =
    # L is the total number of unique points in the simplex grid, which we know a priori
    let L: int = binom(ndiv+dim-1, dim-1)
    result = newTensor[int]([L, dim])
    var x = zeros[int](dim)
    x[dim-1] = ndiv
    for j in 0..dim-1:
        result[0, j] = x[j]
    var h = dim
    for i in 1..L-1:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..dim-1:
            result[i, j] = x[j]
        if val != 1:
            h = dim
    return result

proc simplex_grid_py*(dim: int, ndiv: int): seq[seq[int]] {.exportpy.} = simplex_grid(dim, ndiv).toSeq2D()

proc simplex_grid_fractional*(dim: int,
                              ndiv: int): Tensor[float] =
    result = simplex_grid(dim, ndiv).asType(float)
    result = result.map(x => x / float(ndiv))
    return result

proc simplex_grid_fractional_py*(dim: int, ndiv: int): seq[seq[float]] {.exportpy.} = simplex_grid_fractional(dim, ndiv).toSeq2D()

proc simplex_internal_grid*(dim: int, 
                            ndiv: int): Tensor[int] =

    # L is the total number of unique points inside the simplex grid, which we know a priori
    let L: int = binom(ndiv-1, dim-1)
    result = newTensor[int]([L, dim])
    var x = ones[int](dim)
    x[dim-1] = ndiv+1-dim
    for j in 0..dim-1:
        result[0, j] = x[j]
    var h = dim
    for i in 1..L-1:
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

    result = simplex_internal_grid(dim, ndiv).asType(float)
    result = result.map(x => x / float(ndiv))
    return result

proc simplex_sampling_hed(dim: int,
                          samples: int): Tensor[float] =
    let hypercubesample = randomTensor[float](
        [samples, dim], 
        1.0
        ).map(x => -ln(x))
    let sums = hypercubesample.sum(axis=1)
    result = hypercubesample /. sums

template benchmark(benchmarkName: string, code: untyped) =
    block:
        let t0 = epochTime()
        code
        echo "Task Time: " & $initDuration(microseconds = ((epochTime() - t0)*1e6).int) & "\n"

proc echoHelp*() = echo """

To run nimplex please either (1) provide no arguments and follow the prompts or 
(2) use "-c" or "--config" to provide the configuration per instructions below:

- Provide the 3-letter configuration for task type:
    1. Grid type or uniform random sampling:
        - F: Full grid (including the simplex boundary)
        - I: Internal grid (only points inside the simplex)
        - R: Random uniform sampling using hypercube exponential sampling.
    2. Fractional or Integer positions:
        - F: Fractional grid (points are normalized to fractions of 1)
        - I: Integer grid (points are integers)
    3. Print full result, its shape, or persist in a file:
        - P: Print full grid (presents full result as a table)
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

proc outFunction(config: string, dim: int, ndiv: int, npyName: string, result: Tensor) =
    case config[2]:
        of 'P': echo "Full grid:", result
        of 'N': result.write_npy(npyName)
        else: discard #return nothing, just print the size
    echo "Full sampling shape:", result.shape

proc taskRouterGrid(config: string, dim: int, ndiv: int, npyName: string) =
    case config[0..1]:
        of "FF": outFunction(config, dim, ndiv, npyName, 
                             simplex_grid_fractional(dim, ndiv))
        of "FI": outFunction(config, dim, ndiv, npyName, 
                             simplex_grid(dim, ndiv))
        of "IF": outFunction(config, dim, ndiv, npyName, 
                             simplex_internal_grid_fractional(dim, ndiv))
        of "II": outFunction(config, dim, ndiv, npyName, 
                             simplex_internal_grid(dim, ndiv))
        of "RF": outFunction(config, dim, ndiv, npyName, 
                             simplex_sampling_hed(dim, samples=ndiv))
        else:
            echo "\n--> Invalid configuration in the first 2 config letters."
            quit(1)

proc configValidation(config: string) = 
    assert config.len == 3, "\n--> Invalid configuration lenght. Must be 3 letters."
    assert config[0] in @['F', 'I', 'R'], "\n--> Invalid configuration (in the 1st letter). Must be F, I or R for Full grid, Internal grid, or Random uniform sampling respectively"
    assert config[1] in @['F', 'I'], "\n--> Invalid configuration (in the 2nd letter). Must be F, or I for Fractional positions, or Integer positions respectively"
    if config[0] == 'R':
        assert config[1] == 'F', "\n--> Integer positions not implemented for Random sampling. Must be F for Fractional positions."
    assert config[2] in @['P', 'S', 'N'], "\n--> Invalid configuration (in the 3rd letter). Must be P, S or N for Print full result, Shape, or persist Numpy output respectively"

proc nDivValidation(config: string, nDiv: int, dim: int) = 
    if config[0] == 'I':
        assert ndiv >= dim, "\n--> Invalid number of divisions. Must be greater or equal to the simplex dimension to produce a non-empty internal grid."
    else:
        assert ndiv > 0, "\n--> Invalid number of divisions. Must be a positive integer."

when isMainModule:
    let args = commandLineParams()

    # Interactive
    if args.len == 0:
        echo "Configuration (Full/Internal/Random)(Fractional/Integer)(Print/Shape/Numpysave) - e.g. FFS or R:"
        let config = readLine(stdin)
        configValidation(config)

        echo "Simplex dimensions:"
        let dim = readLine(stdin).parseInt() 

        var nDiv: int
        if config[0]=='R':
            echo "Number of samples:"
            nDiv = readLine(stdin).parseInt()
            assert nDiv > 0, "\n--> Invalid number of samples. Must be a positive integer"
        else:
            echo "N divisions:"
            ndiv = readLine(stdin).parseInt() 
            nDivValidation(config, ndiv, dim)

        var npyName: string = "nimplex_" & config[0..1] & "_" & $dim & "_" & $ndiv & ".npy"
        if config[2] == 'N':
            echo "NumPy Array Filename (skip for default: " & npyName & "):"
            let tempIn = readLine(stdin)
            if tempIn.len > 0:
                npyName = tempIn
            echo "Persisting to NumPy array file:", npyName

        taskRouterGrid(config, dim, ndiv, npyName)

    # Configured
    elif args[0] == "-c" or args[0] == "--config":
        let config = args[1]
        echo "Running with configuration:", args[1..<args.len]
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
            echo "Persisting to NumPy array file:", npyName

        taskRouterGrid(config, dim, ndiv, npyName)

    elif args[0] in @["-h", "--help"]:
        echoHelp()
        quit(0)

    elif args[0] in @["-b", "--benchmark"]:
        benchmark "Simplex Grid Full (dim=9, ndiv=12):":
            discard simplex_grid(9, 12)
        benchmark "Simplex Grid Internal (dim=9, ndiv=12):":
            discard simplex_internal_grid(9, 12)
        benchmark "Simplex Random Sampling (dim=9, samples=1M):":
            discard simplex_sampling_hed(9, 1_000_000)

    # Fallback
    else:
        echoHelp()
        quit(1)