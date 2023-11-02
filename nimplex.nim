# Copyrigth (C) 2023 Adam M. Krajewski

from std/math import binom, ln
import std/sugar
import std/os

import arraymancer
import strutils

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


proc simplex_grid_fractional*(dim: int,
                              ndiv: int): Tensor[float] =

    result = simplex_grid(dim, ndiv).asType(float)
    result = result.map(x => x / float(ndiv))
    return result

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

proc echoHelp*() = echo """

To run the program either (1) provide no arguments and follow the prompts or 
(2) use "-c" or "--config" to provide the configuration.

- For uniform random sampling, provide "R" and:
    -c R [simplex dimension] [number of samples]

- For a grid, provide the 3-letter configuration:
    1. Full or Internal grid:
        - F: Full grid (including the simplex boundary)
        - I: Internal grid (only points inside the simplex)
    2. Fractional or Integer grid:
        - F: Fractional grid (points are normalized to fractions of 1)
        - I: Integer grid (points are integers)
    3. Full or Shape:
        - F: Full grid (present the full result)
        - S: Shape (only the shape / size information)
        - N: Persist to NumPy array file ("nimplex_<configFlags>.npy" or 
             custom path as additonal argument)

    followed by 2 integers of simplex dimensions and number of divisions, like:
        -c FFF [simplex dimension] [number of divisions]
        -c IIF [simplex dimension] [number of divisions]
        -c FIN [simplex dimension] [number of divisions] [path/to/file.npy]
    
"""

proc taskRouter(config: string, dim: int, ndiv: int, npyName: string) =
    let mainConfig = config[0..1]
    assert config[2] == 'F' or config[2] == 'S', "Invalid configuration (in the 3rd letter)"
    case mainConfig:
        of "FF":
            let temp = simplex_grid_fractional(dim, ndiv)
            if config[2] == 'F':
                echo "Full grid:", temp
            echo "Full grid size:", temp.shape[0]
        of "FI":
            let temp = simplex_grid(dim, ndiv)
            if config[2] == 'F':
                echo "Full grid:", temp
            echo "Full grid size:", temp.shape[0]
        of "IF":
            let temp = simplex_internal_grid_fractional(dim, ndiv)
            if config[2] == 'F':
                echo "Full grid:", temp
            echo "Full grid size:", temp.shape[0]
        of "II":
            let temp = simplex_internal_grid(dim, ndiv)
            if config[2] == 'F':
                echo "Full grid:", temp
            echo "Full grid size:", temp.shape[0]
        else:
            echo "Invalid configuration (in the first 2 letters)"
            quit(1)

when isMainModule:
    let args = commandLineParams()

    # Interactive
    if args.len == 0:
        echo "Configuration (Full/Internal/Random)(Fractional/Integer)(Full/Shape) - e.g. FFS or R:"
        let config = readLine(stdin)

        assert config.len == 3 or config=="R", "Invalid configuration"

        echo "Simplex dimensions:"
        let dim = readLine(stdin).parseInt() 

        if config[0]=='R':
            echo "Number of samples:"
            let sampleN = readLine(stdin).parseInt()
            assert sampleN > 0, "Invalid sample number. Must be a positive integer"
            echo simplex_sampling_hed(dim, sampleN)
            quit(0)

        echo "N divisions:"
        let ndiv = readLine(stdin).parseInt() 

        var npyName = "nimplex_" & config[0..1] & "_" & $dim & "_" & $ndiv & ".npy"
        if config[2] == 'N':
            echo "NumPy Array Filename (skip for default: " & npyName & "):"
            let tempIn = readLine(stdin)
            if tempIn.len > 0:
                npyName = tempIn
            echo "Persisting to NumPy array file:", npyName

        taskRouter(config, dim, ndiv, npyName)

    # Configured
    elif args[0] == "-c" or args[0] == "--config":
        let config = args[1]
        assert config.len == 3 or config=="R", "Invalid configuration"

        let dim = args[2].parseInt()
        assert dim > 0, "Invalid dimension"

        if config[0]=='R':
            let sampleN = args[3].parseInt()
            assert sampleN > 0, "Invalid sample number. Must be a positive integer"
            echo simplex_sampling_hed(dim, sampleN)
            quit(0)
        
        let ndiv = args[3].parseInt()

        taskRouter(config, dim, ndiv)

    elif args[0] == "-h" or args[0] == "--help":
        echoHelp()
        quit(0)

    # Fallback
    else:
        echoHelp()
        quit(1)