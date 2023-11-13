# Copyrigth (C) 2023 Adam M. Krajewski
import nimpy

from std/math import binom, ln
import std/sugar

import arraymancer/Tensor 

proc simplex_grid_nim(
    dim: int, 
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

proc simplex_grid*(dim: int, ndiv: int): seq[seq[int]] {.exportpy.} = simplex_grid_nim(dim, ndiv).toSeq2D()

proc simplex_grid_fractional_nim(dim: int,
                              ndiv: int): Tensor[float] =

    result = simplex_grid_nim(dim, ndiv).asType(float)
    result = result.map(x => x / float(ndiv))
    return result

proc simplex_grid_fractional*(dim: int, ndiv: int): seq[seq[float]] {.exportpy.} = simplex_grid_fractional_nim(dim, ndiv).toSeq2D()

proc simplex_internal_grid_nim(dim: int, 
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

proc simplex_internal_grid*(dim: int, ndiv: int): seq[seq[int]] {.exportpy.} = simplex_internal_grid_nim(dim, ndiv).toSeq2D()

proc simplex_internal_grid_fractional_nim(dim: int,
                                       ndiv: int): Tensor[float] =

    result = simplex_internal_grid_nim(dim, ndiv).asType(float)
    result = result.map(x => x / float(ndiv))
    return result

proc simplex_internal_grid_fractional*(dim: int, ndiv: int): seq[seq[float]] {.exportpy.} = simplex_internal_grid_fractional_nim(dim, ndiv).toSeq2D()

proc simplex_sampling_hed_nim(dim: int,
                          samples: int): Tensor[float] =
    let hypercubesample = randomTensor[float](
        [samples, dim], 
        1.0
        ).map(x => -ln(x))
    let sums = hypercubesample.sum(axis=1)
    result = hypercubesample /. sums

proc simplex_sampling_hed*(dim: int, samples: int): seq[seq[float]] {.exportpy.} = simplex_sampling_hed_nim(dim, samples).toSeq2D()

proc simplex_graph2_nim(
    dim: int, 
    ndiv: int): (Tensor[int], seq[seq[int]]) =

    # L is the total number of unique points in the simplex grid, which we know a priori
    let L: int = binom(ndiv+dim-1, dim-1)
    var result1 = newTensor[int]([L, dim])
    var neighbors = newSeq[seq[int]](L)
    var x = zeros[int](dim)

    func neighborsLink(i:int, x:Tensor): seq[int] =
        if x[0] != 0:
            result.add(i-1)
        if x[1] != 0:
            result.add(i+1)
        return result

    x[dim-1] = ndiv
    for j in 0..dim-1:
        result1[0, j] = x[j]
    var h = dim

    neighbors[0] = neighborsLink(0, x)
    echo x.toSeq1D(), neighbors[0], "\n"

    for i in 1..L-1:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..dim-1:
            result1[i, j] = x[j]
        echo x.toSeq1D()
        neighbors[i] = neighborsLink(i, x)
        echo neighbors[i]
        echo "\n"
        if val != 1:
            h = dim
    return (result1, neighbors)

proc simplex_graph2*(dim: int, ndiv: int): (seq[seq[int]], seq[seq[int]]) {.exportpy.} = 
    let graph = simplex_graph2_nim(dim, ndiv)
    return (graph[0].toSeq2D(), graph[1])

proc simplex_graph3_nim(
    dim: int, 
    ndiv: int): (Tensor[int], seq[seq[int]]) =

    # L is the total number of unique points in the simplex grid, which we know a priori
    let L: int = binom(ndiv+dim-1, dim-1)
    var result1 = newTensor[int]([L, dim])
    var neighbors = newSeq[seq[int]](L)
    var x = zeros[int](dim)

    func neighborsLink(i:int, x:Tensor): seq[int] =
        let jump0 = 1
        let jump1 = binom(1+12-x[0], 1)

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

    x[dim-1] = ndiv
    for j in 0..dim-1:
        result1[0, j] = x[j]
    var h = dim

    neighbors[0] = neighborsLink(0, x)

    for i in 1..L-1:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..dim-1:
            result1[i, j] = x[j]
        neighbors[i] = neighborsLink(i, x)
        if val != 1:
            h = dim
    return (result1, neighbors)

proc simplex_graph3*(dim: int, ndiv: int): (seq[seq[int]], seq[seq[int]]) {.exportpy.} = 
    let graph = simplex_graph3_nim(dim, ndiv)
    return (graph[0].toSeq2D(), graph[1])

proc simplex_graph3_fractional_nim(dim: int,
                                  ndiv: int): (Tensor[float], seq[seq[int]]) =
    let tempResult = simplex_graph3_nim(dim, ndiv)
    return (tempResult[0].asType(float).map(x => x / float(ndiv)), tempResult[1])

proc simplex_graph3_fractional*(dim: int, ndiv: int): (seq[seq[float]], seq[seq[int]]) {.exportpy.} = 
    let graph = simplex_graph3_fractional_nim(dim, ndiv)
    return (graph[0].toSeq2D(), graph[1])

proc simplex_graph4_nim(
    dim: int, 
    ndiv: int): (Tensor[int], seq[seq[int]]) =

    # L is the total number of unique points in the simplex grid, which we know a priori
    let L: int = binom(ndiv+dim-1, dim-1)
    var result1 = newTensor[int]([L, dim])
    var neighbors = newSeq[seq[int]](L)
    var x = zeros[int](dim)

    func neighborsLink(i:int, x:Tensor): seq[int] =
        let jump0 = 1  #binom(12-x[0]-x[1]-x[2], 0)
        let jump1 = binom(1+12-x[0]-x[1], 1)
        let jump2 = binom(2+12-x[0], 2)
        if x[0] != 0:
            result.add(i - jump2)                 # quaternary
            result.add(i - jump2 - jump1)         # quaternary
            result.add(i - jump2 - jump1 - jump0) # quaternary

        if x[1] != 0:
            result.add(i - jump1)                 # ternary
            result.add(i - jump1 - jump0)         # ternary
            result.add(i + jump2 - jump1 - x[1])  # quaternary

        if x[2] != 0:
            result.add(i - jump0)                 # binary
            result.add(i + jump1 - jump0)         # ternary
            result.add(i + jump2 - jump0 - x[1])  # quaternary

        if x[3] != 0:
            result.add(i + jump0)                 # binary
            result.add(i + jump1)                 # ternary
            result.add(i + jump2 - x[1])          # quaternary  

        return result

    x[dim-1] = ndiv
    for j in 0..dim-1:
        result1[0, j] = x[j]
    var h = dim

    neighbors[0] = neighborsLink(0, x)
    #echo x.toSeq1D(), neighbors[0], "\n"

    for i in 1..L-1:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..dim-1:
            result1[i, j] = x[j]
        neighbors[i] = neighborsLink(i, x)
        #echo x.toSeq1D(), neighbors[i],  "\n"
        if val != 1:
            h = dim
    return (result1, neighbors)

proc simplex_graph4*(dim: int, ndiv: int): (seq[seq[int]], seq[seq[int]]) {.exportpy.} = 
    let graph = simplex_graph4_nim(dim, ndiv)
    return (graph[0].toSeq2D(), graph[1])

proc simplex_graph5_nim(
    dim: int, 
    ndiv: int): (Tensor[int], seq[seq[int]]) =

    # L is the total number of unique points in the simplex grid, which we know a priori
    let L: int = binom(ndiv+dim-1, dim-1)
    var result1 = newTensor[int]([L, dim])
    var neighbors = newSeq[seq[int]](L)
    var x = zeros[int](dim)

    func neighborsLink(i:int, x:Tensor): seq[int] =
        let jump0 = 1  #binom(12-x[0]-x[1]-x[2]-x[3], 0)
        let jump1 = binom(1+12-x[0]-x[1]-x[2], 1)
        let jump2 = binom(2+12-x[0]-x[1], 2)

        if x[1] != 0:
            result.add(i - jump2)                 # quaternary
            result.add(i - jump2 - jump1)         # quaternary
            result.add(i - jump2 - jump1 - jump0) # quaternary

        if x[2] != 0:
            result.add(i - jump1)                 # ternary
            result.add(i - jump1 - jump0)         # ternary
            result.add(i + jump2 - jump1 - x[2])  # quaternary

        if x[3] != 0:
            result.add(i - jump0)                 # binary
            result.add(i + jump1 - jump0)         # ternary
            result.add(i + jump2 - jump0 - x[2])  # quaternary

        if x[4] != 0:
            result.add(i + jump0)                 # binary
            result.add(i + jump1)                 # ternary
            result.add(i + jump2 - x[2])          # quaternary  

        return result

    x[dim-1] = ndiv
    for j in 0..dim-1:
        result1[0, j] = x[j]
    var h = dim

    neighbors[0] = neighborsLink(0, x)
    #echo x.toSeq1D(), neighbors[0], "\n"

    for i in 1..L-1:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..dim-1:
            result1[i, j] = x[j]
        neighbors[i] = neighborsLink(i, x)
        #echo x.toSeq1D(), neighbors[i],  "\n"
        if val != 1:
            h = dim
    return (result1, neighbors)

proc simplex_graph5*(dim: int, ndiv: int): (seq[seq[int]], seq[seq[int]]) {.exportpy.} = 
    let graph = simplex_graph5_nim(dim, ndiv)
    return (graph[0].toSeq2D(), graph[1])