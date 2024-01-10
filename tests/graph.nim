import std/unittest
import ../nimplex
import arraymancer/Tensor
import std/sequtils
import std/sugar
import std/sets

# SMALL GRAPHS
echo "*** SMALL GRAPHS ***"

suite "small simplex integer 2-component (binary) graph":
    let (nodes, neighbors) = nimplex.simplex_graph(2, 5)
    let neighborsNumber: seq[int] = neighbors.map(n => n.len)
        
    echo "Nodes:\n", nodes
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 2
    test "correct number of nodes/vertices":
        check nodes.shape[0] == 6
    test "correct maximum number of neighbors":
        check neighborsNumber.max == 2*(2-1)
    test "correct minimum number of neighbors":
        check neighborsNumber.min == (2-1)
    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == 
            @[@[0, 5], @[1, 4], @[2, 3], @[3, 2], @[4, 1], @[5, 0]]
    test "correct neighbors list for each node/vertex":
        check neighbors == 
            @[@[1], @[0, 2], @[1, 3], @[2, 4], @[3, 5], @[4]]

suite "small simplex integer 3-component (ternary) graph":
    let (nodes, neighbors) = nimplex.simplex_graph(3, 3)
    let neighborsNumber: seq[int] = neighbors.map(n => n.len)
        
    echo "Nodes:\n", nodes
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 3
    test "correct number of nodes/vertices":
        check nodes.shape[0] == 10
    test "correct maximum number of neighbors":
        check neighborsNumber.max == 3*(3-1)
    test "correct minimum number of neighbors":
        check neighborsNumber.min == (3-1)
    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == @[
            @[0, 0, 3], @[0, 1, 2], @[0, 2, 1], @[0, 3, 0], 
            @[1, 0, 2], @[1, 1, 1], @[1, 2, 0], 
            @[2, 0, 1], @[2, 1, 0], 
            @[3, 0, 0]]
    test "correct neighbors list for each node/vertex":
        check neighbors == @[
            @[1, 4], @[0, 2, 4, 5], @[1, 3, 5, 6], @[2, 6], 
            @[1, 0, 5, 7], @[2, 1, 4, 6, 7, 8], @[3, 2, 5, 8], 
            @[5, 4, 8, 9], @[6, 5, 7, 9], 
            @[8, 7]]

suite "small simplex integer 4-component (quaternary) graph":
    let (nodes, neighbors) = nimplex.simplex_graph(4, 4)
    let neighborsNumber: seq[int] = neighbors.map(n => n.len)

    echo "Nodes:\n", nodes
    echo "Neighbors:\n", neighbors

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 4
    test "correct number of nodes/vertices":
        check nodes.shape[0] == 35
    test "correct maximum number of neighbors":
        check neighborsNumber.max == 4*(4-1)
    test "correct minimum number of neighbors":
        check neighborsNumber.min == (4-1)
    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == @[
            @[0, 0, 0, 4], @[0, 0, 1, 3], @[0, 0, 2, 2], @[0, 0, 3, 1], @[0, 0, 4, 0], 
            @[0, 1, 0, 3], @[0, 1, 1, 2], @[0, 1, 2, 1], @[0, 1, 3, 0], 
            @[0, 2, 0, 2], @[0, 2, 1, 1], @[0, 2, 2, 0], 
            @[0, 3, 0, 1], @[0, 3, 1, 0], 
            @[0, 4, 0, 0], 
            @[1, 0, 0, 3], @[1, 0, 1, 2], @[1, 0, 2, 1], @[1, 0, 3, 0], 
            @[1, 1, 0, 2], @[1, 1, 1, 1], @[1, 1, 2, 0], 
            @[1, 2, 0, 1], @[1, 2, 1, 0], 
            @[1, 3, 0, 0], 
            @[2, 0, 0, 2], @[2, 0, 1, 1], @[2, 0, 2, 0], 
            @[2, 1, 0, 1], @[2, 1, 1, 0], 
            @[2, 2, 0, 0], 
            @[3, 0, 0, 1], @[3, 0, 1, 0], 
            @[3, 1, 0, 0], 
            @[4, 0, 0, 0]]
    test "correct neighbors list for each node/vertex":
        check neighbors == @[
            @[1, 5, 15], @[0, 2, 5, 6, 15, 16], @[1, 3, 6, 7, 16, 17], @[2, 4, 7, 8, 17, 18], @[3, 8, 18], 
            @[1, 0, 6, 9, 15, 19], @[2, 1, 5, 7, 9, 10, 16, 19, 20], @[3, 2, 6, 8, 10, 11, 17, 20, 21], @[4, 3, 7, 11, 18, 21], 
            @[6, 5, 10, 12, 19, 22], @[7, 6, 9, 11, 12, 13, 20, 22, 23], @[8, 7, 10, 13, 21, 23], 
            @[10, 9, 13, 14, 22, 24], @[11, 10, 12, 14, 23, 24], 
            @[13, 12, 24], 
            @[5, 1, 0, 16, 19, 25], @[6, 2, 1, 15, 17, 19, 20, 25, 26], @[7, 3, 2, 16, 18, 20, 21, 26, 27], @[8, 4, 3, 17, 21, 27], 
            @[9, 6, 5, 16, 15, 20, 22, 25, 28], @[10, 7, 6, 17, 16, 19, 21, 22, 23, 26, 28, 29], @[11, 8, 7, 18, 17, 20, 23, 27, 29], 
            @[12, 10, 9, 20, 19, 23, 24, 28, 30], @[13, 11, 10, 21, 20, 22, 24, 29, 30],
            @[14, 13, 12, 23, 22, 30], 
            @[19, 16, 15, 26, 28, 31], @[20, 17, 16, 25, 27, 28, 29, 31, 32], @[21, 18, 17, 26, 29, 32], 
            @[22, 20, 19, 26, 25, 29, 30, 31, 33], @[23, 21, 20, 27, 26, 28, 30, 32, 33], 
            @[24, 23, 22, 29, 28, 33],
            @[28, 26, 25, 32, 33, 34], @[29, 27, 26, 31, 33, 34], 
            @[30, 29, 28, 32, 31, 34], 
            @[33, 32, 31]]
        
suite "small simplex fractional 4-component (quaternary) graph":
    # Only testing the node positions, since the neighbors are the same as in the integer case suite above
    let (nodes, _) = nimplex.simplex_graph_fractional(4, 4)
    echo "Nodes:\n", nodes

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 4
    test "correct number of nodes/vertices":
        check nodes.shape[0] == 35
    test "correct node/vertex positions in the simplex":
        check nodes.toSeq2D() == @[
            @[0.0, 0.0, 0.0, 1.0], @[0.0, 0.0, 0.25, 0.75], @[0.0, 0.0, 0.5, 0.5], @[0.0, 0.0, 0.75, 0.25], @[0.0, 0.0, 1.0, 0.0], 
            @[0.0, 0.25, 0.0, 0.75], @[0.0, 0.25, 0.25, 0.5], @[0.0, 0.25, 0.5, 0.25], @[0.0, 0.25, 0.75, 0.0], 
            @[0.0, 0.5, 0.0, 0.5], @[0.0, 0.5, 0.25, 0.25], @[0.0, 0.5, 0.5, 0.0], 
            @[0.0, 0.75, 0.0, 0.25], @[0.0, 0.75, 0.25, 0.0], 
            @[0.0, 1.0, 0.0, 0.0], 
            @[0.25, 0.0, 0.0, 0.75], @[0.25, 0.0, 0.25, 0.5], @[0.25, 0.0, 0.5, 0.25], @[0.25, 0.0, 0.75, 0.0], 
            @[0.25, 0.25, 0.0, 0.5], @[0.25, 0.25, 0.25, 0.25], @[0.25, 0.25, 0.5, 0.0], 
            @[0.25, 0.5, 0.0, 0.25], @[0.25, 0.5, 0.25, 0.0],
            @[0.25, 0.75, 0.0, 0.0],
            @[0.5, 0.0, 0.0, 0.5], @[0.5, 0.0, 0.25, 0.25], @[0.5, 0.0, 0.5, 0.0],
            @[0.5, 0.25, 0.0, 0.25], @[0.5, 0.25, 0.25, 0.0],
            @[0.5, 0.5, 0.0, 0.0],
            @[0.75, 0.0, 0.0, 0.25], @[0.75, 0.0, 0.25, 0.0],
            @[0.75, 0.25, 0.0, 0.0],
            @[1.0, 0.0, 0.0, 0.0]]


# LARGE GRAPHS 
echo "*** LARGE GRAPHS ***"

suite "large simplex integer 3-component (ternary) graph":
    let 
        ndiv = 48
        (nodes, neighbors) = nimplex.simplex_graph(3, ndiv)
        neighborsNumber: seq[int] = neighbors.map(n => n.len)

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 3
    test "correct number of nodes/vertices":
        check nodes.shape[0] == 1225
    test "correct maximum number of neighbors":
        check neighborsNumber.max == 3*(3-1)
    test "correct minimum number of neighbors":
        check neighborsNumber.min == (3-1)
    test "correct positions in the simplex of a cherry-picked node at index 0":
        check nodes[0, _].toSeq2D()[0] == @[0, 0, ndiv]
    test "correct positions in the simplex of a cherry-picked node at index 123":
        check nodes[123, _].toSeq2D()[0] == @[2, 26, 20]
    test "correct neighbors list for a cherry-picked node at index 0":
        check neighbors[0] == @[1, ndiv+1]
    test "correct neighbors list for a cherry-picked node at index 123":
        check neighbors[123] == @[76, 75, 122, 124, 169, 170]
    
suite "3C special case method agreement with general case (20k nodes)":
    let
        ndiv = 200
        (nodes1, neighbors1) = nimplex.simplex_graph(3, ndiv)
        (nodes2, neighbors2) = nimplex.simplex_graph_3C(ndiv)

    test "correct dimensionality of nodes/vertices for both methods":
        check nodes1.shape[1] == 3
        check nodes2.shape[1] == nodes1.shape[1]
    test "correct number of nodes/vertices for both methods":
        check nodes1.shape[0] == 20301
        check nodes2.shape[0] == nodes1.shape[0]
    test "both methods produce the same nodes/vertices":
        check nodes1 == nodes2
    test "both methods produce the same neighbors (in any order))":
        for i in 0..<nodes1.shape[0]:
            check neighbors1[i].toHashSet() == neighbors2[i].toHashSet()

suite "very large simplex fractional 12-component graph (1M+ nodes  / 9M+ edges)":
    let 
        ndiv = 12
        (nodes, neighbors) = nimplex.simplex_graph(12, ndiv)
        neighborsNumbers: seq[int] = neighbors.map(n => n.len)
        edgesCount = neighborsNumbers.foldl(a+b)

    test "correct dimensionality of nodes/vertices":
        check nodes.shape[1] == 12

    test "correct number of nodes/vertices":
        check nodes.shape[0] == 1_352_078

    test "correct number of neighbors (graph edges)":
        check edgesCount == 93_117_024

    test "correct maximum number of neighbors":
        check neighborsNumbers.max == 12*(12-1)

    test "correct minimum number of neighbors":
        check neighborsNumbers.min == (12-1)

    test "correct positions in the simplex of a cherry-picked node at index 0":
        check nodes[0, _].toSeq2D()[0] == @[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ndiv]

    test "correct positions in the simplex of a cherry-picked node at index 123_456":
        check nodes[123_456, _].toSeq2D()[0] == @[0, 0, 0, 6, 0, 1, 0, 0, 1, 1, 1, 2]

    test "correct neighbors list for a cherry-picked node at index 0":
        check neighbors[0] == @[1, 13, 91, 455, 1820, 6188, 18564, 50388, 125970, 293930, 646646]
        
    test "correct neighbors list for a cherry-picked node at index 123_456":
        check neighbors[123_456] == @[
            121740, 120816, 120564, 120438, 120382, 120367, 120363, 120362, 123204, 123078, 123022, 123007, 
            123003, 123002, 123441, 123437, 123436, 123452, 123451, 123455, 123457, 123459, 123460, 123466, 
            123469, 123470, 123491, 123501, 123504, 123505, 123561, 123571, 123574, 123575, 123687, 123697, 
            123700, 123701, 123918, 124149, 124159, 124162, 124163, 124710, 124941, 124951, 124954, 124955, 
            199038, 200292, 200523, 200533, 200536, 200537, 366998, 368252, 368483, 368493, 368496, 368497, 
            719714, 720968, 721199, 721209, 721212, 721213]