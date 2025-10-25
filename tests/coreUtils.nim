import std/unittest
import ../nimplex
import std/times
import std/tables
import std/strutils
import std/math

# SMALL GRIDS
let t0 = cpuTime()
echo "***** CORE UTILS *****"

proc twoHotEdgePrintout(edgeKey: seq[int], components: seq[seq[float]]): string =
    let alphabets = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for pos in 0..<edgeKey.len:
        var compStr: string = ""
        if edgeKey[pos] == 1:
            for i in 0..<components[pos].len:
                let v = components[pos][i]
                if v != 0.0:
                    compStr.add(alphabets[i] & $int(v*100) & " ")
            result.add(compStr & " <-> ")
    result.removeSuffix(" <-> ")


suite "Pure components scaling sanity test":
    const components: seq[seq[float]] = @[
            @[1.0, 0.0, 0.0], # Pure Fe
            @[0.0, 1.0, 0.0], # Pure C
            @[0.0, 0.0, 1.0]  # Pure Mo
    ]

    let testTable: Table[seq[int], float] = nimplex.elemental_space_scaling_table(components)
    for key, value in testTable.pairs:
        test "Scaling for pure component " & $key & " -> " & twoHotEdgePrintout(key, components) & " is approximately 1.0 as expected":
            check almostEqual(value, 1.0)

suite "Simple case of scaling (30-60-90 triangle, but in simplex space)":
    const 
        components: seq[seq[float]] = @[
            @[1.0, 0.0, 0.0], # A
            @[0.5, 0.5, 0.0], # AB
            @[0.0, 0.0, 1.0]  # C
        ]
        referenceTable: Table[seq[int], float] = {
            @[1, 1, 0]: 0.5,       # A-AB edge
            @[0, 1, 1]: 1.0,         # AB-C edge
            @[1, 0, 1]: 1.0        # A-C edge
        }.toTable

    let testTable: Table[seq[int], float] = nimplex.elemental_space_scaling_table(components)
    for key, value in testTable.pairs:
        test "Scaling for " & $key & " -> " & twoHotEdgePrintout(key, components) & " matches reference (" & $referenceTable[key] & ")":
            check almostEqual(value, referenceTable[key])

suite "Small real elemental space scaling table test":
    const 
        components: seq[seq[float]] = @[
            @[0.94, 0.05, 0.01], # Fe95 C5 Mo1
            @[3.0, 1.0, 0.0],    # Fe3C
            @[0.2, 0.0, 0.8]     # Fe20 Mo80
        ]
        referenceTable: Table[seq[int], float] = {
            @[0, 1, 1]: 0.8, 
            @[1, 1, 0]: 0.2, 
            @[1, 0, 1]: 0.79
        }.toTable

    let testTable = nimplex.elemental_space_scaling_table(components)
    for key, value in testTable.pairs:
        test "Scaling for " & $key & " -> " & twoHotEdgePrintout(key, components) & " matches reference (" & $referenceTable[key] & ")":
            check almostEqual(value, referenceTable[key])

echo "***** CORE UTILS BENCHMARK RESULTS *****"
let t1 = cpuTime()
echo "Time:\n" & $initDuration(milliseconds = ((t1 - t0)*1e3).int) & "\n"
echo "---------------------------------\n"
