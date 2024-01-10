import std/unittest
import ../nimplex
import arraymancer/Tensor

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