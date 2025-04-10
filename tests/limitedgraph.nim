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

