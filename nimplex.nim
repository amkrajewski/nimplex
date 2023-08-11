# Copyrigth (C) 2023 Adam M. Krajewski

from std/math import binom
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

echo "Simplex dimensions:"
let dim = readLine(stdin).parseInt() 

echo "N divisions:"
let ndiv = readLine(stdin).parseInt() 

let temp = simplex_internal_grid(dim, ndiv)
echo "Internal grid:", temp
echo "Internal grid size:", temp.shape[0]
