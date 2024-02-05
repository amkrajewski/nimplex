import std/unittest
import ../src/nimplex
import arraymancer/Tensor
import std/times

let t0 = cpuTime()
# SMALL GRIDS

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
        # Nice Finite Decimal Fractions we can check for exact equality directly
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
        # Infinite Decimal Fractions, so we check for approximate equality
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

suite "verify attainable simplex grid given by pure components is equivalent to the simplex grid itself (fractional)":
    let result0 = simplex_grid_fractional(3,6)
    let result1 = result0.attainable2elemental(@[@[1.0,0.0,0.0],@[0.0,1.0,0.0],@[0.0,0.0,1.0]])
    test "matching shape":
        check result0.shape == result1.shape
    test "matching values":
        for i in 0..27:
            for j in 0..2:
                check abs(result0[i, j] - result1[i, j]) < 0.0001

suite "small simplex attainable grids (binary with ndiv=6 in ternary from [1,1,1] to [1,0,0])":
    let result = simplex_grid_fractional(2,6).attainable2elemental(@[@[1.0,0.0,0.0],@[1.0,1.0,1.0]])
    echo "Shape:", result.shape
    echo "Result:", result.toSeq2D()
    test "correct dimensionality":
        check result.shape[1] == 3
    test "grid has correct number of nodes/vertices":
        check result.shape[0] == 7
    test "correct values":
        let reference = @[
            @[6/18, 6/18, 6/18],
            @[8/18, 5/18, 5/18],
            @[10/18, 4/18, 4/18],
            @[12/18, 3/18, 3/18],
            @[14/18, 2/18, 2/18],
            @[16/18, 1/18, 1/18],
            @[18/18, 0/18, 0/18]]
        for i in 0..6:
            for j in 0..2:
                check abs(result[i, j] - reference[i][j]) < 0.0001

let t1 = cpuTime()
# LARGE GRIDS

suite "large simplex full integer grids":
    let result = nimplex.simplex_grid(7, 24)
    let resultCherryPick = result[124_492, _].toSeq2D()[0]
    echo "Shape:  ", result.shape
    echo "Result #124_492:  ", resultCherryPick
    test "correct dimensionality":
        check result.shape[1] == 7
    test "grid has correct number of nodes/vertices":
        check result.shape[0] == 593_775
    test "correct value at index 124_492":
        check resultCherryPick == @[1, 0, 2, 3, 9, 8, 1]
    
suite "small simplex internal integer grids":
    let result = nimplex.simplex_internal_grid(7, 24)
    let resultCherryPick = result[44_234, _].toSeq2D()[0]
    echo "Shape:  ", result.shape
    echo "Result #44_234:  ", resultCherryPick
    test "correct dimensionality":
        check result.shape[1] == 7
    test "grid has correct number of vertices":
        check result.shape[0] == 100_947
    test "correct value at index 44_234":
        check resultCherryPick == @[2, 7, 3, 2, 1, 4, 5]

suite "large simplex full fractional grids":
    let result = nimplex.simplex_grid_fractional(7, 24)
    let resultCherryPick = result[124_492, _].toSeq2D()[0]
    echo "Shape:  ", result.shape
    echo "Result #124_492:  ", resultCherryPick
    test "correct dimensionality":
        check result.shape[1] == 7
    test "grid has correct number of nodes/vertices":
        check result.shape[0] == 593_775
    test "correct value at index 124_492":
        # Infinite Decimal Fractions, so we check for approximate equality
        let reference = @[0.041666, 0.0, 0.083333, 0.125, 0.375, 0.333333, 0.041666]
        for i in 0..6:
            check abs(resultCherryPick[i] - reference[i]) < 0.0001

suite "large simplex internal fractional grids":
    let result = nimplex.simplex_internal_grid_fractional(7, 24)
    let resultCherryPick = result[44_234, _].toSeq2D()[0]
    echo "Shape:  ", result.shape
    echo "Result #44_234:  ", resultCherryPick
    test "correct dimensionality":
        check result.shape[1] == 7
    test "grid has correct number of vertices":
        check result.shape[0] == 100_947
    test "correct value at index 44_234":
        # Infinite Decimal Fractions, so we check for approximate equality
        let reference = @[0.0833333, 0.291666, 0.125, 0.0833333, 0.0416666, 0.166666, 0.208333]
        for i in 0..6:
            check abs(resultCherryPick[i] - reference[i]) < 0.0001

let t2 = cpuTime()
echo "Small Grids:\n" & $initDuration(microseconds = ((t1 - t0)*1e6).int) & "\n" 
echo "Large Grids:\n" & $initDuration(milliseconds = ((t2 - t1)*1e3).int) & "\n"
