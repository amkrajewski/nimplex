import std/unittest
import ../nimplex
import arraymancer/Tensor
import std/sequtils
import std/sugar
import std/sets
import std/times
import std/math

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

suite "medium (12-divisions) simplex integer 3-component (ternary) graph with limits with limits of [[1, 7], [1, 7], [1, 7]] for a symmetrical limit around center":
    let 
        nDiv:int = 12
        limit:seq[seq[int]] = @[@[1, nDiv-5], @[1, nDiv-5], @[1, nDiv-5]]
        (nodes, neighbors) = nimplex.simplex_graph_limited(3, nDiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
    
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 3

    test "correct number of nodes/vertices":
        let tempNodes = nimplex.simplex_grid(3, nDiv)
        # Here, let's manually count the number of nodes in the grid that are within the limits, as visual inspection is less reliable.
        var count = 0
        for i in 0..<tempNodes.shape[0]:
            var withinLimits = true
            for j in 0..<tempNodes.shape[1]:
                if tempNodes[i, j] < limit[j][0] or tempNodes[i, j] > limit[j][1]:
                    withinLimits = false
                    break
            if withinLimits:
                count += 1
        check nodes.shape[0] == count

    test "correct number of neighbors (graph edges)":
        check edgesCount == 180

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 3*(3-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == 3

    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == 
            @[
                @[1, 4, 7], 
                @[1, 5, 6], @[1, 6, 5], @[1, 7, 4], @[2, 3, 7], @[2, 4, 6], @[2, 5, 5], @[2, 6, 4], @[2, 7, 3], @[3, 2, 7], @[3, 3, 6], @[3, 4, 5], @[3, 5, 4], @[3, 6, 3], 
                @[3, 7, 2], @[4, 1, 7], @[4, 2, 6], @[4, 3, 5], @[4, 4, 4], @[4, 5, 3], @[4, 6, 2], @[4, 7, 1], @[5, 1, 6], @[5, 2, 5], @[5, 3, 4], @[5, 4, 3], @[5, 5, 2], 
                @[5, 6, 1], @[6, 1, 5], @[6, 2, 4], @[6, 3, 3], @[6, 4, 2], @[6, 5, 1], @[7, 1, 4], @[7, 2, 3], @[7, 3, 2], @[7, 4, 1]
            ]

    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[
                @[1, 4, 5], @[0, 2, 5, 6], @[1, 3, 6, 7], @[2, 7, 8], @[0, 5, 9, 10], @[1, 0, 4, 6, 10, 11], @[2, 1, 5, 7, 11, 12], @[3, 2, 6, 8, 12, 13], @[3, 7, 13, 14], 
                @[4, 10, 15, 16], @[5, 4, 9, 11, 16, 17], @[6, 5, 10, 12, 17, 18], @[7, 6, 11, 13, 18, 19], @[8, 7, 12, 14, 19, 20], @[8, 13, 20, 21], @[9, 16, 22], 
                @[10, 9, 15, 17, 22, 23], @[11, 10, 16, 18, 23, 24], @[12, 11, 17, 19, 24, 25], @[13, 12, 18, 20, 25, 26], @[14, 13, 19, 21, 26, 27], @[14, 20, 27], 
                @[16, 15, 23, 28], @[17, 16, 22, 24, 28, 29], @[18, 17, 23, 25, 29, 30], @[19, 18, 24, 26, 30, 31], @[20, 19, 25, 27, 31, 32], @[21, 20, 26, 32], 
                @[23, 22, 29, 33], @[24, 23, 28, 30, 33, 34], @[25, 24, 29, 31, 34, 35], @[26, 25, 30, 32, 35, 36], @[27, 26, 31, 36], @[29, 28, 34], @[30, 29, 33, 35], 
                @[31, 30, 34, 36], @[32, 31, 35]
            ]

suite "medium (12-divisions) simplex integer 3-component (ternary) graph with limits with limits of [[0, 7], [0, 7], [0, 7]] for an asymmetrical limit around center":
    let 
        nDiv:int = 12
        limit:seq[seq[int]] = @[@[0, 7], @[0, 7], @[0, 7]]
        (nodes, neighbors) = nimplex.simplex_graph_limited(3, nDiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
    
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 3

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 46

    test "correct number of neighbors (graph edges)":
        check edgesCount == 228

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 3*(3-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == 3

    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == 
            @[
                @[0, 5, 7], @[0, 6, 6], @[0, 7, 5], @[1, 4, 7], @[1, 5, 6], @[1, 6, 5], @[1, 7, 4], @[2, 3, 7], @[2, 4, 6], @[2, 5, 5], @[2, 6, 4], @[2, 7, 3], @[3, 2, 7], 
                @[3, 3, 6], @[3, 4, 5], @[3, 5, 4], @[3, 6, 3], @[3, 7, 2], @[4, 1, 7], @[4, 2, 6], @[4, 3, 5], @[4, 4, 4], @[4, 5, 3], @[4, 6, 2], @[4, 7, 1], @[5, 0, 7], 
                @[5, 1, 6], @[5, 2, 5], @[5, 3, 4], @[5, 4, 3], @[5, 5, 2], @[5, 6, 1], @[5, 7, 0], @[6, 0, 6], @[6, 1, 5], @[6, 2, 4], @[6, 3, 3], @[6, 4, 2], @[6, 5, 1], 
                @[6, 6, 0], @[7, 0, 5], @[7, 1, 4], @[7, 2, 3], @[7, 3, 2], @[7, 4, 1], @[7, 5, 0]
            ]

    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[
                @[1, 3, 4], @[0, 2, 4, 5], @[1, 5, 6], @[0, 4, 7, 8], @[1, 0, 3, 5, 8, 9], @[2, 1, 4, 6, 9, 10], @[2, 5, 10, 11], @[3, 8, 12, 13], @[4, 3, 7, 9, 13, 14], 
                @[5, 4, 8, 10, 14, 15], @[6, 5, 9, 11, 15, 16], @[6, 10, 16, 17], @[7, 13, 18, 19], @[8, 7, 12, 14, 19, 20], @[9, 8, 13, 15, 20, 21], @[10, 9, 14, 16, 21, 22], 
                @[11, 10, 15, 17, 22, 23], @[11, 16, 23, 24], @[12, 19, 25, 26], @[13, 12, 18, 20, 26, 27], @[14, 13, 19, 21, 27, 28], @[15, 14, 20, 22, 28, 29], 
                @[16, 15, 21, 23, 29, 30], @[17, 16, 22, 24, 30, 31], @[17, 23, 31, 32], @[18, 26, 33], @[19, 18, 25, 27, 33, 34], @[20, 19, 26, 28, 34, 35], 
                @[21, 20, 27, 29, 35, 36], @[22, 21, 28, 30, 36, 37], @[23, 22, 29, 31, 37, 38], @[24, 23, 30, 32, 38, 39], @[24, 31, 39], @[26, 25, 34, 40], 
                @[27, 26, 33, 35, 40, 41], @[28, 27, 34, 36, 41, 42], @[29, 28, 35, 37, 42, 43], @[30, 29, 36, 38, 43, 44], @[31, 30, 37, 39, 44, 45], @[32, 31, 38, 45], 
                @[34, 33, 41], @[35, 34, 40, 42], @[36, 35, 41, 43], @[37, 36, 42, 44], @[38, 37, 43, 45], @[39, 38, 44]
            ]

suite "medium (12-divisions) simplex integer 3-component (ternary) graph with limits with limits of [[1, 8], [1, 8], [1, 8]] for an asymmetrical limit around center":
    let 
        nDiv:int = 12
        limit:seq[seq[int]] = @[@[1, 8], @[1, 8], @[1, 8]]
        (nodes, neighbors) = nimplex.simplex_graph_limited(3, nDiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
    
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 3

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 46

    test "correct number of neighbors (graph edges)":
        check edgesCount == 228

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 3*(3-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == 3

    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == 
            @[
                @[1, 3, 8], @[1, 4, 7], @[1, 5, 6], @[1, 6, 5], @[1, 7, 4], @[1, 8, 3], @[2, 2, 8], @[2, 3, 7], @[2, 4, 6], @[2, 5, 5], @[2, 6, 4], @[2, 7, 3], @[2, 8, 2], @[3, 1, 8], 
                @[3, 2, 7], @[3, 3, 6], @[3, 4, 5], @[3, 5, 4], @[3, 6, 3], @[3, 7, 2], @[3, 8, 1], @[4, 1, 7], @[4, 2, 6], @[4, 3, 5], @[4, 4, 4], @[4, 5, 3], @[4, 6, 2], @[4, 7, 1], 
                @[5, 1, 6], @[5, 2, 5], @[5, 3, 4], @[5, 4, 3], @[5, 5, 2], @[5, 6, 1], @[6, 1, 5], @[6, 2, 4], @[6, 3, 3], @[6, 4, 2], @[6, 5, 1], @[7, 1, 4], @[7, 2, 3], @[7, 3, 2], 
                @[7, 4, 1], @[8, 1, 3], @[8, 2, 2], @[8, 3, 1]
            ]

    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[
                @[1, 6, 7], @[0, 2, 7, 8], @[1, 3, 8, 9], @[2, 4, 9, 10], @[3, 5, 10, 11], @[4, 11, 12], @[0, 7, 13, 14], @[1, 0, 6, 8, 14, 15], @[2, 1, 7, 9, 15, 16], 
                @[3, 2, 8, 10, 16, 17], @[4, 3, 9, 11, 17, 18], @[5, 4, 10, 12, 18, 19], @[5, 11, 19, 20], @[6, 14, 21], @[7, 6, 13, 15, 21, 22], @[8, 7, 14, 16, 22, 23], 
                @[9, 8, 15, 17, 23, 24], @[10, 9, 16, 18, 24, 25], @[11, 10, 17, 19, 25, 26], @[12, 11, 18, 20, 26, 27], @[12, 19, 27], @[14, 13, 22, 28], @[15, 14, 21, 23, 28, 29], 
                @[16, 15, 22, 24, 29, 30], @[17, 16, 23, 25, 30, 31], @[18, 17, 24, 26, 31, 32], @[19, 18, 25, 27, 32, 33], @[20, 19, 26, 33], @[22, 21, 29, 34], 
                @[23, 22, 28, 30, 34, 35], @[24, 23, 29, 31, 35, 36], @[25, 24, 30, 32, 36, 37], @[26, 25, 31, 33, 37, 38], @[27, 26, 32, 38], @[29, 28, 35, 39], 
                @[30, 29, 34, 36, 39, 40], @[31, 30, 35, 37, 40, 41], @[32, 31, 36, 38, 41, 42], @[33, 32, 37, 42], @[35, 34, 40, 43], @[36, 35, 39, 41, 43, 44],
                @[37, 36, 40, 42, 44, 45], @[38, 37, 41, 45], @[40, 39, 44], @[41, 40, 43, 45], @[42, 41, 44]
            ]


suite "medium (12-divisions) simplex integer 4-component (quaternary) graph with limits with limits of [[0, 7], [0, 7], [0, 7], [0, 7]] for an asymmetrical limit around center":
    let 
        nDiv:int = 12
        limit:seq[seq[int]] = @[@[0, 7], @[0, 7], @[0, 7], @[0, 7]]
        (nodes, neighbors) = nimplex.simplex_graph_limited(4, nDiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
        
    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 4

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 315

    test "correct number of neighbors (graph edges)":
        check edgesCount == 3048

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 4*(4-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == 5

    test "correct positions in the simplex of a cherry-picked node at index 0":
        check nodes[0, _].toSeq2D()[0] == @[0, 0, 5, 7]

    test "correct positions in the simplex of a random-picked node at index 123":
        check nodes[123, _].toSeq2D()[0] == @[2, 4, 3, 3]

    test "correct neighbors list for a cherry-picked node at index 0":
        check neighbors[0] == @[1, 3, 4, 46, 47]

    test "correct neighbors list for a random-picked node at index 123":
        check neighbors[123] == @[79, 72, 71, 116, 115, 122, 124, 129, 130, 166, 172, 173]

suite "*fractional* small (10-divisions) simplex integer 2-component (binary) graph with limits with limits of [[0, 1], [0.04, 0.666]]":
    let 
        nDiv: int = 10
        limit: seq[seq[float]] = @[@[0.0, 1.0], @[0.049, 0.666]]
        (nodes, neighbors) = nimplex.simplex_graph_limited_fractional(2, nDiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
    
    echo "Nodes:\n", nodes
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 2

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 8

    test "correct number of neighbors (graph edges)":
        check edgesCount == 14

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 2*(2-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == (2-1)

    test "correct node/vertex positions in the simplex":
        let refSeq: seq[seq[float]] = 
            @[
                @[0.3, 0.7], 
                @[0.4, 0.6], 
                @[0.5, 0.5], 
                @[0.6, 0.4], 
                @[0.7, 0.3], 
                @[0.8, 0.2], 
                @[0.9, 0.1], 
                @[1.0, 0.0]
            ]
        check nodes.toSeq2D() == refSeq
            

    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[@[1], @[0, 2], @[1, 3], @[2, 4], @[3, 5], @[4, 6], @[5, 7], @[6]]

suite "*fractional* medium (12-divisions) simplex integer 4-component (quaternary) graph with limits with limits of [0.0123, 0.58] for an asymmetrical limit around center, which should agree with integer limit of [0, 7]":
    let 
        nDiv:int = 12
        limit:seq[seq[float]] = @[@[0.0123, 0.58], @[0.0123, 0.58], @[0.0123, 0.58], @[0.0123, 0.58]]
        (nodes, neighbors) = nimplex.simplex_graph_limited_fractional(4, nDiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)
        
    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 4

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 315

    test "correct number of neighbors (graph edges)":
        check edgesCount == 3048

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 4*(4-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == 5

    test "correct positions in the simplex of a cherry-picked node at index 0":
        for pair in zip(nodes[0, _].toSeq2D()[0], @[0.0, 0.0, 0.4167, 0.5833]):
            check round(pair[0], 3) == round(pair[1], 3)

    test "correct positions in the simplex of a random-picked node at index 123":
        for pair in zip(nodes[123, _].toSeq2D()[0], @[0.1667, 0.3333, 0.25, 0.25]):
            check round(pair[0], 3) == round(pair[1], 3)

    test "correct neighbors list for a cherry-picked node at index 0":
        check neighbors[0] == @[1, 3, 4, 46, 47]

    test "correct neighbors list for a random-picked node at index 123":
        check neighbors[123] == @[79, 72, 71, 116, 115, 122, 124, 129, 130, 166, 172, 173]

let t1 = cpuTime()

suite "very large simplex 12-component graph (1M nodes  / 67M edges) with limits of [0, 4] (i.e. up to 33%) in each component":
    let 
        ndiv: int = 12
        limit: seq[seq[int]] = collect: 
            for i in 0..<ndiv: 
                @[0, 4]
        (nodes, neighbors) = nimplex.simplex_graph_limited(12, ndiv, limit)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 12

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 975_338

    test "correct number of neighbors (graph edges)":
        check edgesCount == 68_704_416

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 12*(12-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == 27

    test "correct positions in the simplex of a cherry-picked node at index 0":
        check nodes[0, _].toSeq2D()[0] == @[0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4]

    test "correct positions in the simplex of a random-picked node at index 123_456":
        check nodes[123_456, _].toSeq2D()[0] == @[0, 0, 2, 0, 1, 0, 0, 1, 1, 3, 3, 1]

    test "correct neighbors list for a cherry-picked node at index 0":
        check neighbors[0] == 
            @[1, 2, 3, 35, 36, 37, 320, 321, 322, 1751, 1752, 1753, 7140, 7141, 7142, 23940, 23941, 23942, 69675, 69676, 69677, 182005, 182006, 182007, 436348, 436349, 436350]

    test "correct neighbors list for a random-picked node at index 123_456":
        check neighbors[123_456] == 
            @[
                92897, 80632, 77343, 76218, 75905, 75835, 75821, 75817, 75816, 120167, 119042, 118729, 118659, 118645, 118641, 118640, 123386, 123372, 123368, 123367, 123442, 
                123438, 123437, 123452, 123451, 123455, 123457, 123460, 123461, 123472, 123476, 123477, 123526, 123542, 123546, 123547, 123741, 123811, 123827, 123831, 123832, 
                124622, 124692, 124708, 124712, 124713, 126911, 126981, 126997, 127001, 127002, 131801, 135256, 135326, 135342, 135346, 135347, 150872, 154327, 154397, 154413, 
                154417, 154418, 235786, 263202, 266657, 266727, 266743, 266747, 266748, 490129, 517545, 521000, 521070, 521086, 521090, 521091
            ]

let t2 = cpuTime()

echo "\n***** LIMITED GRAPH BENCHMARK RESULTS *****\n"
echo "Small Graphs:\n" & $initDuration(microseconds = ((t1 - t0)*1e6).int) & "\n"
echo "Large Graphs:\n" & $initDuration(milliseconds = ((t2 - t1)*1e3).int) & "\n"
echo "-------------------------------------\n"
