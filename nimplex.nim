# Copyrigth (C) 2023 Adam M. Krajewski

from std/math import binom, ln
import std/sugar
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

when isMainModule:
    echo "Configuration (Full/Internal)(Fractional/Integer)(Full/Shape) - e.g. FFS:"
    let config = readLine(stdin)
    
    echo "Simplex dimensions:"
    let dim = readLine(stdin).parseInt() 

    echo "N divisions:"
    let ndiv = readLine(stdin).parseInt() 

    let mainConfig = config[0..1]

    case mainConfig
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
        echo "Invalid configuration"
        quit(1)
