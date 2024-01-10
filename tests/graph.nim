import std/unittest
import ../nimplex
import arraymancer/Tensor

suite "small simplex integer 2-component (binary) graph":
    let (nodes, neighbors) = nimplex.simplex_graph(2, 5)
    var neighborsNumber: seq[int]
    for i in 0..5:
        neighborsNumber.add(neighbors[i].len)
        
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
    var neighborsNumber: seq[int]
    for i in 0..5:
        neighborsNumber.add(neighbors[i].len)
        
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