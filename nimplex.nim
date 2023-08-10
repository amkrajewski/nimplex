# Copyrigth (C) 2023 Adam M. Krajewski

import std/math
import arraymancer
import strutils

proc simplex_grid*(m: int, 
                   n: int): Tensor[int] =
    let L: int = binom(n+m-1, m-1)
    result = newTensor[int]([L, m])
    var x = zeros[int](m)
    x[m-1] = n
    for j in 0..m-1:
        result[0, j] = x[j]
    var h = m
    for i in 1..L-1:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[m-1] = val - 1
        x[h-1] += 1
        for j in 0..m-1:
            result[i, j] = x[j]
        if val != 1:
            h = m
    return result

echo "Simplex dimensions:"
let m = readLine(stdin).parseInt() 

echo "N divisions:"
let n = readLine(stdin).parseInt() 

echo simplex_grid(m, n)