import std/unittest
import ../nimplex
import arraymancer/Tensor
import std/sequtils
import std/sugar

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
        
    