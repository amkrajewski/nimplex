import std/tables
import std/unittest
import std/times
import std/strformat
import std/random
import std/sequtils
import arraymancer/Tensor

import ../nimplex
import ../utils/stitching

echo "***** Stitching (Graph Complex Construction) Tests *****"

suite "Set up two 3C simplex grids that happen to share a common 2C subspace and verify correct stitching":
    let t0 = cpuTime()
    # Set up A B C in A B C D
    let grid1 = nimplex.simplex_grid_fractional(3, 5).attainable2elemental(@[
            @[1.0,0.0,0.0,0.0],
            @[0.0,1.0,0.0,0.0],
            @[0.0,0.0,1.0,0.0]]
            )

    # Set up D C A in A B C D
    let grid2 = nimplex.simplex_grid_fractional(3, 5).attainable2elemental(@[
            @[0.0,0.0,0.0,1.0],
            @[0.0,0.0,1.0,0.0],
            @[1.0,0.0,0.0,0.0]]
            )

    test "Points generated and are in different spaces":
        check grid1 != grid2

    let stitch1Table = findStitchingPoints(3, 5, components = @["A", "B", "C"])
    let stitch2Table = findStitchingPoints(3, 5, components = @["D", "C", "A"])

    test "Stitching tables generated and of correct size (6 ternaries, 6 binaries, 3 unaries = 15)":
        check stitch1Table.len == 15
        check stitch2Table.len == 15

    let stitch1 = stitch1Table["C-A"]
    let stitch2 = stitch2Table["C-A"]

    test "Binary stitching points were extracted from table and are of correct size (ndiv+1 = 6)":
        check stitch1.len == 6
        check stitch2.len == 6

    test "Binary stitching points are differently indexed (because they are in different spaces)":
        check stitch1 != stitch2

    test "Stitching points match expected reference values":
        check stitch1 == @[0, 6, 11, 15, 18, 20]
        check stitch2 == @[5, 4, 3, 2, 1, 0]

    for i in 0..<6:
        let p1 = grid1[stitch1[i], _].squeeze().toSeq1D()
        let p2 = grid2[stitch2[i], _].squeeze().toSeq1D()
        test fmt"Stitching point {i} ({p1}) in grid1 matches stitching point {i} in grid2 ({p2})":
            check p1 == p2

    echo fmt"Finished small test in {round(1000*(cpuTime()-t0), 4)} ms."

suite "Enhance the previous suite by (1) assigning higher order elemental space with chemical element names (Ti Al V Cr Ni) \nand (2) setting up the elemental space to be composed of alloys rather than pure elements":
    let comps = @["Ti", "Al", "V", "Cr", "Ni"]
    let t0 = cpuTime()
    # A - (Ti0.8 Al0.15 V0.05) TiAlloy
    # B - V0.25 Cr0.75 - CrV
    # C - V0.75 Ni0.25 - VNi
    # D - Al0.9 Cr0.1 - AlCr

    # Set up A B C 
    let grid1 = nimplex.simplex_grid_fractional(3, 5).attainable2elemental(@[
        #   Ti    Al     V    Cr    Ni
        @[ 0.8, 0.15, 0.05,  0.0,  0.0], # A
        @[ 0.0,  0.0, 0.25, 0.75,  0.0], # B
        @[ 0.0,  0.0, 0.75, 0.0,  0.25]] # C
        )

    # Set up D C A
    let grid2 = nimplex.simplex_grid_fractional(3, 5).attainable2elemental(@[
        #   Ti    Al     V    Cr    Ni
        @[ 0.0,  0.9,  0.0,  0.1,  0.0], # D
        @[ 0.0,  0.0, 0.75,  0.0, 0.25], # C
        @[ 0.8, 0.15, 0.05,  0.0,  0.0]] # A
        )

    test "Points generated and are in different spaces":
        check grid1 != grid2

    let stitch1Table = findStitchingPoints(3, 5, components = @["TiAlloy", "CrV", "VNi"])
    let stitch2Table = findStitchingPoints(3, 5, components = @["AlCr", "VNi", "TiAlloy"])

    test "Stitching tables generated and of correct size (6 ternaries, 6 binaries, 3 unaries = 15)":
        check stitch1Table.len == 15
        check stitch2Table.len == 15

    let stitch1 = stitch1Table["TiAlloy-VNi"]
    let stitch2 = stitch2Table["TiAlloy-VNi"]

    test "Binary stitching points were extracted from table and are of correct size (ndiv+1 = 5)":
        check stitch1.len == 6
        check stitch2.len == 6

    test "Binary stitching points are differently indexed (because they are in different spaces)":
        check stitch1 != stitch2

    test "Stitching points match expected reference values":
        check stitch1 == @[20, 18, 15, 11, 6, 0]
        check stitch2 == @[0, 1, 2, 3, 4, 5]

    func formula(fracs: seq[float], components: seq[string]): string = 
        for (frac, comp) in zip(fracs, components):
            if frac > 0:
                result.add(fmt"{comp}{frac:.2f} ")

    for i in 0..5:
        let p1 = grid1[stitch1[i], _].squeeze().toSeq1D().mapIt(round(it, 3))
        let p2 = grid2[stitch2[i], _].squeeze().toSeq1D().mapIt(round(it, 3))
        test fmt"Stitching point {i} ({p1:<30}) in grid1 matches stitching point {i} in grid2 ({p2:<30}) -> {formula(p1, comps)}":
            check p1 == p2

    echo fmt"Finished medium test in {round(1000*(cpuTime()-t0), 4)} ms."


suite "Set up dissimilar 6C and 7C simplex grids in 9C elemental space that happen to share a common 4C subspace and verify correct stitching":
    test "Don't try to visualize this, it's too big! :)":
        check true

    let t0 = cpuTime()

    # Let's say A D H I are common between the two spaces

    # Inside A B C D E F G H I set up A G D H I B (shared + G B)
    let grid1: Tensor[float] = nimplex.simplex_grid(6, 9).attainable2elemental(@[
            # A  B  C  D  E  F  G  H  I
            @[1, 0, 0, 0, 0, 0, 0, 0, 0],
            @[0, 0, 0, 0, 0, 0, 1, 0, 0],
            @[0, 0, 0, 1, 0, 0, 0, 0, 0],
            @[0, 0, 0, 0, 0, 0, 0, 1, 0],
            @[0, 0, 0, 0, 0, 0, 0, 0, 1],
            @[0, 1, 0, 0, 0, 0, 0, 0, 0]]
            )

    # Inside A B C D E F G H I set up H I E F A C D (shared + E F C)
    let grid2: Tensor[float] = nimplex.simplex_grid(7, 9).attainable2elemental(@[
            # A  B  C  D  E  F  G  H  I
            @[0, 0, 0, 0, 0, 0, 0, 1, 0],
            @[0, 0, 0, 0, 0, 0, 0, 0, 1],
            @[0, 0, 0, 0, 1, 0, 0, 0, 0],
            @[0, 0, 0, 0, 0, 1, 0, 0, 0],
            @[1, 0, 0, 0, 0, 0, 0, 0, 0],
            @[0, 0, 1, 0, 0, 0, 0, 0, 0],
            @[0, 0, 0, 1, 0, 0, 0, 0, 0]]
            )

    test fmt"Generated spaces are of different sizes ({grid1.shape} != {grid2.shape})":
        check grid1.shape != grid2.shape

    test "Generated spaces are different":
        check grid1 != grid2

    let stitch1Table = findStitchingPoints(6, 9, maxDim=4, components = @["A", "G", "D", "H", "I", "B"])
    let stitch2Table = findStitchingPoints(7, 9, maxDim=4, components = @["H", "I", "E", "F", "A", "C", "D"])

    test "Stitching table 1 is of correct size 516 for <=4C in 6C space":
        check stitch1Table.len == 516

    test "Stitching table 2 is of correct size 1099 for <=4C in 7C space":
        check stitch2Table.len == 1099

    let stitch1 = stitch1Table["A-D-H-I"]
    let stitch2 = stitch2Table["A-D-H-I"]

    test "Binary stitching points were extracted from table and are of correct size: binom(ndiv+dim-1, dim-1) or in this case 220":
        check stitch1.len == 220
        check stitch2.len == 220

    test "Binary stitching points are differently indexed (because they are in different spaces)":
        check stitch1 != stitch2

    test "The first 25 stitching points match expected reference values":
        check stitch1[0..<25] == @[2001, 1999, 1998, 1997, 1990, 1989, 1988, 1986, 1985, 1983, 1965, 1964, 1963, 1961, 1960, 1958, 1955, 1954, 1952, 1949, 1910, 1909, 1908, 1906, 1905]
        check stitch2[0..<25] == @[54, 52, 2046, 759, 49, 2044, 757, 3324, 2532, 1245, 45, 2041, 754, 3322, 2530, 1243, 4108, 3646, 2854, 1567, 40, 2037, 750, 3319, 2527]

    # Check 15 random points in the stitching tables
    let pointsToVerify: seq[int] = newSeq[int](15).mapIt(rand(220))

    for i in 0..<15:
        let p1 = grid1[stitch1[pointsToVerify[i]], _].squeeze().toSeq1D()
        let p2 = grid2[stitch2[pointsToVerify[i]], _].squeeze().toSeq1D()
        test fmt"Stitching Grid 1 point {pointsToVerify[i]:<4} ({p1}) matches Stitching Grid 2 point {pointsToVerify[i]:<4} in grid2 ({p2})":
            check p1 == p2

    echo fmt"Finished large test in {round(1000*(cpuTime()-t0), 1)} ms."

