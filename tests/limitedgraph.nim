import std/unittest
import ../nimplex
import arraymancer/Tensor
import std/sequtils
import std/sugar
import std/sets
import std/times

# SMALL GRAPHS
let t0 = cpuTime()
echo "***** SMALL TESTS *****"

suite "small (5-divisions) simplex integer 2-component (binary) graph with limits with limits of [[0, 5], [0, 5]] so *no limit*":
    let 
        nDiv:int = 5
        limit:seq[seq[int]] = @[@[0, nDiv], @[0, nDiv]]
        (nodes, neighbors) = nimplex.simplex_graph_limited(2, 5, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
    
    echo "Nodes:\n", nodes
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 2

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 6

    test "correct number of neighbors (graph edges)":
        check edgesCount == 10

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 2*(2-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == (2-1)

    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == 
            @[@[0, 5], @[1, 4], @[2, 3], @[3, 2], @[4, 1], @[5, 0]]

    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[@[1], @[0, 2], @[1, 3], @[2, 4], @[3, 5], @[4]]

suite "small (5-divisions) simplex integer 2-component (binary) graph with limits with limits of [[0, 5], [0, 3]] so *limit on single (second) component*":
    let 
        nDiv:int = 5
        limit:seq[seq[int]] = @[@[0, nDiv], @[0, 3]]
        (nodes, neighbors) = nimplex.simplex_graph_limited(2, 5, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
    
    echo "Nodes:\n", nodes
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 2

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 4

    test "correct number of neighbors (graph edges)":
        check edgesCount == 6

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 2*(2-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == 1

    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == 
            @[@[2, 3], @[3, 2], @[4, 1], @[5, 0]]

    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[@[1], @[0, 2], @[1, 3], @[2]]


suite "small (5-divisions) simplex integer 3-component (ternary) graph with limits with limits of [[0, 5], [0, 5], [0, 3]] so limit on third component":
    let 
        nDiv:int = 5
        limit:seq[seq[int]] = @[@[0, nDiv], @[0, nDiv], @[0, 3]]
        (nodes, neighbors) = nimplex.simplex_graph_limited(3, nDiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
    
    echo "Nodes:\n", nodes
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 3

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 18

    test "correct number of neighbors (graph edges)":
        check edgesCount == 76

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 3*(3-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == (3-1)

    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == 
            @[
                @[0, 2, 3], 
                @[0, 3, 2], 
                @[0, 4, 1], 
                @[0, 5, 0], 
                @[1, 1, 3], 
                @[1, 2, 2], 
                @[1, 3, 1], 
                @[1, 4, 0], 
                @[2, 0, 3], 
                @[2, 1, 2], 
                @[2, 2, 1], 
                @[2, 3, 0], 
                @[3, 0, 2], 
                @[3, 1, 1], 
                @[3, 2, 0], 
                @[4, 0, 1], 
                @[4, 1, 0], 
                @[5, 0, 0]
            ]

    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[
                @[1, 4, 5], 
                @[0, 2, 5, 6], 
                @[1, 3, 6, 7], 
                @[2, 7], 
                @[0, 5, 8, 9], 
                @[1, 0, 4, 6, 9, 10], 
                @[2, 1, 5, 7, 10, 11], 
                @[3, 2, 6, 11], 
                @[4, 9, 12], 
                @[5, 4, 8, 10, 12, 13], 
                @[6, 5, 9, 11, 13, 14], 
                @[7, 6, 10, 14], 
                @[9, 8, 13, 15], 
                @[10, 9, 12, 14, 15, 16], 
                @[11, 10, 13, 16], 
                @[13, 12, 16, 17], 
                @[14, 13, 15, 17], 
                @[16, 15]
            ]

suite "medium (12-divisions) simplex integer 3-component (ternary) graph with limits with limits of [[0, 12], [0, 5], [2, 10]] so limit two components":
    let 
        nDiv:int = 12
        limit:seq[seq[int]] = @[@[0, nDiv], @[0, 5], @[2, 10]]
        (nodes, neighbors) = nimplex.simplex_graph_limited(3, nDiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
    
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 3

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 48

    test "correct number of neighbors (graph edges)":
        check edgesCount == 236

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 3*(3-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == (3-1)

    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == 
            @[
                @[0, 2, 10], @[0, 3, 9], @[0, 4, 8], @[0, 5, 7], @[1, 1, 10], @[1, 2, 9], @[1, 3, 8], @[1, 4, 7], @[1, 5, 6], @[2, 0, 10], @[2, 1, 9], @[2, 2, 8], @[2, 3, 7], @[2, 4, 6], 
                @[2, 5, 5], @[3, 0, 9], @[3, 1, 8], @[3, 2, 7], @[3, 3, 6], @[3, 4, 5], @[3, 5, 4], @[4, 0, 8], @[4, 1, 7], @[4, 2, 6], @[4, 3, 5], @[4, 4, 4], @[4, 5, 3], @[5, 0, 7], 
                @[5, 1, 6], @[5, 2, 5], @[5, 3, 4], @[5, 4, 3], @[5, 5, 2], @[6, 0, 6], @[6, 1, 5], @[6, 2, 4], @[6, 3, 3], @[6, 4, 2], @[7, 0, 5], @[7, 1, 4], @[7, 2, 3], @[7, 3, 2], 
                @[8, 0, 4], @[8, 1, 3], @[8, 2, 2], @[9, 0, 3], @[9, 1, 2], @[10, 0, 2]
            ]

    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[
                @[1, 4, 5], @[0, 2, 5, 6], @[1, 3, 6, 7], @[2, 7, 8], @[0, 5, 9, 10], @[1, 0, 4, 6, 10, 11], @[2, 1, 5, 7, 11, 12], @[3, 2, 6, 8, 12, 13], @[3, 7, 13, 14], @[4, 10, 15], 
                @[5, 4, 9, 11, 15, 16], @[6, 5, 10, 12, 16, 17], @[7, 6, 11, 13, 17, 18], @[8, 7, 12, 14, 18, 19], @[8, 13, 19, 20], @[10, 9, 16, 21], @[11, 10, 15, 17, 21, 22], 
                @[12, 11, 16, 18, 22, 23], @[13, 12, 17, 19, 23, 24], @[14, 13, 18, 20, 24, 25], @[14, 19, 25, 26], @[16, 15, 22, 27], @[17, 16, 21, 23, 27, 28], @[18, 17, 22, 24, 28, 29], 
                @[19, 18, 23, 25, 29, 30], @[20, 19, 24, 26, 30, 31], @[20, 25, 31, 32], @[22, 21, 28, 33], @[23, 22, 27, 29, 33, 34], @[24, 23, 28, 30, 34, 35], @[25, 24, 29, 31, 35, 36], 
                @[26, 25, 30, 32, 36, 37], @[26, 31, 37], @[28, 27, 34, 38], @[29, 28, 33, 35, 38, 39], @[30, 29, 34, 36, 39, 40], @[31, 30, 35, 37, 40, 41], @[32, 31, 36, 41], 
                @[34, 33, 39, 42], @[35, 34, 38, 40, 42, 43], @[36, 35, 39, 41, 43, 44], @[37, 36, 40, 44], @[39, 38, 43, 45], @[40, 39, 42, 44, 45, 46], @[41, 40, 43, 46], 
                @[43, 42, 46, 47], @[44, 43, 45, 47], @[46, 45]
            ]

