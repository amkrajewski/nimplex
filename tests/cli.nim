import std/[os, osproc]
import std/unittest
import std/strutils
import std/sequtils
import std/sugar
import arraymancer/Tensor
import arraymancer/io
import std/hashes
import ../nimplex

suite "test if correct output is given when nimplex is run in command line with some selected configurations":
    test "check if compiled nimplex is present in the current working directory":
        require fileExists("nimplex")

    test "generate small integer simplex grid (FIS 3 3) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c FIS 3 3")
            outputLines = output.splitLines
            reference = @["Running with configuration:@[\"FIS\", \"3\", \"3\"]", "Full shape:[10, 3]"]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i] == reference[i]

    test "generate large integer simplex grid (FIS 9 12) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c FIS 9 12")
            outputLines = output.splitLines
            reference = @["Running with configuration:@[\"FIS\", \"9\", \"12\"]", "Full shape:[125970, 9]"]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i] == reference[i]

    test "generate large internal integer simplex grid (IIS 7 12) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c IIS 7 12")
            outputLines = output.splitLines
            reference = @["Running with configuration:@[\"IIS\", \"7\", \"12\"]", "Full shape:[462, 7]"]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i] == reference[i]

    test "generate small integer simplex grid (FIP 3 3) and print values to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c FIP 3 3")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration:@[\"FIP\", \"3\", \"3\"]", 
                "Full Output:Tensor[system.int] of shape \"[10, 3]\" on backend \"Cpu\"", 
                "|0      0     3|", 
                "|0      1     2|", 
                "|0      2     1|", 
                "|0      3     0|", 
                "|1      0     2|", 
                "|1      1     1|", 
                "|1      2     0|", 
                "|2      0     1|", 
                "|2      1     0|", 
                "|3      0     0|", 
                "Full shape:[10, 3]"]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i] == reference[i]

    test "generate small integer internal simplex grid (IIP 3 5) and print values to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c IIP 3 5")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration:@[\"IIP\", \"3\", \"5\"]", 
                "Full Output:Tensor[system.int] of shape \"[6, 3]\" on backend \"Cpu\"", 
                "|1      1     3|", 
                "|1      2     2|", 
                "|1      3     1|", 
                "|2      1     2|", 
                "|2      2     1|", 
                "|3      1     1|", 
                "Full shape:[6, 3]"]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i] == reference[i]

    test "(significant if previous failes) values themself generated for small integer internal simplex grid (IIP 3 5) corectness":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c IIP 3 5")
            outputLines = output.splitLines
            reference = @[
                @[1, 1, 3],
                @[1, 2, 2],
                @[1, 3, 1],
                @[2, 1, 2],
                @[2, 2, 1],
                @[3, 1, 1]]
        check exitCode == 0
        var outputGrid = newSeq[seq[int]](outputLines.len-3)
        for i in 0..<reference.len:
            for v in outputLines[i+2].replace("|","").splitWhitespace:
                outputGrid[i].add(parseInt(v))

        for i in 0..<reference.len:
            check outputGrid[i] == reference[i]

    test "(independent of numerical precision of a given machine) fractional internal simplex grid (IFP 3 6) correctness":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c IFP 3 6")
            outputLines = output.splitLines
            reference = @[
                @[0.166667, 0.166667, 0.666667],
                @[0.166667, 0.333333, 0.5],
                @[0.166667, 0.5, 0.333333],
                @[0.166667, 0.666667, 0.166667],
                @[0.333333, 0.166667, 0.5],
                @[0.333333, 0.333333, 0.333333],
                @[0.333333, 0.5, 0.166667],
                @[0.5, 0.166667, 0.333333],
                @[0.5, 0.333333, 0.166667],
                @[0.666667, 0.166667, 0.166667]]

        check exitCode == 0
        var outputGrid = newSeq[seq[float]](outputLines.len-3)
        for i in 0..<reference.len:
            for v in outputLines[i+2].replace("|","").splitWhitespace:
                outputGrid[i].add(parseFloat(v))

        for i in 0..<reference.len:
            check abs(outputGrid[i][0] - reference[i][0]) < 0.0001

suite "Test NumPy exports corectness":

    test "generate auto-named a medium fractional internal simplex grid (IFP 7 11) and export it to NumPy (nimplex_IF_7_11.npy)":
        let (output, exitCode) = execCmdEx("./nimplex -c IFN 7 11")
        check exitCode == 0
        check fileExists("nimplex_IF_7_11.npy")

    test "generate a medium fractional internal simplex grid (IFP 7 11) named testExport.npy and export it to NumPy":
        let (output, exitCode) = execCmdEx("./nimplex -c IFN 7 11 testExport.npy")
        check exitCode == 0
        check fileExists("testExport.npy")

    test "Verify that hashes of the two NumPy exports are exactly the same binary":
        let 
            hash1 = readFile("nimplex_IF_7_11.npy").hash.toHex
            hash2 = readFile("testExport.npy").hash.toHex
        check hash1 == hash2

    test "Verify that the exported NumPy files match the direct result from nimplex (tested for corectness in grid.nim)":
        let 
            result = nimplex.simplex_internal_grid_fractional(7, 11)
            loadedNumpy = read_npy[float]("testExport.npy")

        check result.shape == loadedNumpy.shape
        check result == loadedNumpy


    