import std/tables
import std/unittest
import std/times
import std/strformat
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

    echo fmt"Finished test in {cpuTime()-t0} seconds."
    