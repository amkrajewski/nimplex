import std/unittest
import ../nimplex
import arraymancer/Tensor
import math

suite "small simplex full integer grids":
    let result = nimplex.simplex_grid(2, 5)
    echo "Shape:", result.shape
    echo "Result:", result.toSeq2D()
    test "correct dimensionality":
        check result.shape[1] == 2
    test "grid has correct number of nodes/vertices":
        check result.shape[0] == 6
    test "correct values":
        check result.toSeq2D() == 
            @[@[0, 5], @[1, 4], @[2, 3], @[3, 2], @[4, 1], @[5, 0]]
        

suite "small simplex internal integer grids":
    let result = nimplex.simplex_internal_grid(2, 7)
    echo "Shape:", result.shape
    echo "Result:", result.toSeq2D()
    test "correct dimensionality":
        check result.shape[1] == 2
    test "grid has correct number of vertices":
        check result.shape[0] == 6
    test "correct values":
        check result.toSeq2D() ==
            @[@[1, 6], @[2, 5], @[3, 4], @[4, 3], @[5, 2], @[6, 1]]

suite "small simplex full fractional grids":
    let result = nimplex.simplex_grid_fractional(2, 5)
    echo "Shape:", result.shape
    echo "Result:", result.toSeq2D()
    test "correct dimensionality":
        check result.shape[1] == 2
    test "grid has correct number of nodes/vertices":
        check result.shape[0] == 6
    test "correct values":
        check result.toSeq2D() == 
            @[@[0.0, 1.0], @[0.2, 0.8], @[0.4, 0.6], @[0.6, 0.4], @[0.8, 0.2], @[1.0, 0.0]]

suite "small simplex internal fractional grids":
    let result = nimplex.simplex_internal_grid_fractional(2, 7)
    echo "Shape:", result.shape
    echo "Result:", result.toSeq2D()
    test "correct dimensionality":
        check result.shape[1] == 2
    test "grid has correct number of vertices":
        check result.shape[0] == 6
    test "correct values":
        let reference = @[
            @[0.142857, 0.857142], 
            @[0.285714, 0.714285], 
            @[0.428571, 0.571428], 
            @[0.571428, 0.428571], 
            @[0.714285, 0.285714], 
            @[0.857142, 0.142857]]
        for i in 0..5:
            for j in 0..1:
                check abs(result[i, j] - reference[i][j]) < 0.0001