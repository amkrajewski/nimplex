# Copyrigth (C) 2023 Adam M. Krajewski

import std/math
import std/sugar
import arraymancer
import strutils

proc simplex_grid*(dim: int, 
                   ndiv: int): Tensor[int] =

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


echo "Simplex dimensions:"
let dim = readLine(stdin).parseInt() 

echo "N divisions:"
let ndiv = readLine(stdin).parseInt() 

echo simplex_grid(dim, ndiv)