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

