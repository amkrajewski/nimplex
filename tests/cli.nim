import std/[os, osproc]
import system
import std/unittest
import std/strutils
import std/sequtils
import arraymancer/Tensor
import arraymancer/io
import std/hashes
import ../nimplex
import std/sets
import std/times
import std/re

let t0 = cpuTime()

proc stripAnsiEscapeCodes(s: string): string =
  let regex = re(r"\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]")
  return s.replace(regex, "")

echo "***** CLI Tests *****"

suite "test if correct grid output is given when nimplex is run in command line with some selected configurations":
    test "check if compiled nimplex is present in the current working directory":
        # For Unix systems, the executable is named nimplex, but for Windows it is nimplex.exe
        echo "Detected host OS: " & hostOS
        if hostOS == "windows":
            require fileExists("nimplex.exe")
        else:
            require fileExists("nimplex")

    test "generate small integer simplex grid (FIS 3 3) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c FIS 3 3")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"FIS\", \"3\", \"3\"]", 
                "Shape: [10, 3]"]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate large integer simplex grid (FIS 9 12) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c FIS 9 12")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"FIS\", \"9\", \"12\"]", 
                "Shape: [125970, 9]"]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate large internal integer simplex grid (IIS 7 12) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c IIS 7 12")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"IIS\", \"7\", \"12\"]",
                "Shape: [462, 7]"]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate small integer simplex grid (FIP 3 3) and print values to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c FIP 3 3")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"FIP\", \"3\", \"3\"]", 
                "Full Output: Tensor[system.int] of shape \"[10, 3]\" on backend \"Cpu\"", 
                "|0      0     3|", 
                "|0      1     2|", 
                "|0      2     1|", 
                "|0      3     0|", 
                "|1      0     2|", 
                "|1      1     1|", 
                "|1      2     0|", 
                "|2      0     1|", 
                "|2      1     0|", 
                "|3      0     0|"
                ]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate small integer internal simplex grid (IIP 3 5) and print values to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c IIP 3 5")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"IIP\", \"3\", \"5\"]", 
                "Full Output: Tensor[system.int] of shape \"[6, 3]\" on backend \"Cpu\"", 
                "|1      1     3|", 
                "|1      2     2|", 
                "|1      3     1|", 
                "|2      1     2|", 
                "|2      2     1|", 
                "|3      1     1|"
                ]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

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
            for j in 0..<3:
                check abs(outputGrid[i][j] - reference[i][j]) < 0.0001

suite "test if correct graph output is given when nimplex is run in command line with some selected configurations":

    test "check if compiled nimplex is present in the current working directory":
        # For Unix systems, the executable is named nimplex, but for Windows it is nimplex.exe
        if hostOS == "windows":
            require fileExists("nimplex.exe")
        else:
            require fileExists("nimplex")

    test "generate small integer simplex graph (GIS 3 3) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c GIS 3 3")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"GIS\", \"3\", \"3\"]", 
                "Nodes Shape: [10, 3]",
                "Edges Count: 36"
                ]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate small limited integer simplex graph (LIS 3 3 [[0,3],[1,2],[0,3]]) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c LIS 3 3 \"[[0,3],[1,2],[0,3]]\"")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"LIS\", \"3\", \"3\", \"[[0,3],[1,2],[0,3]]\"]", 
                "Nodes Shape: [5, 3]",
                "Edges Count: 14"
                ]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate small limited fractional simplex graph (LFS 3 3 [[0,1],[0.333,0.666],[0,1]]) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c LFS 3 3 \"[[0,1],[0.333,0.666],[0,1]]\"")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"LFS\", \"3\", \"3\", \"[[0,1],[0.333,0.666],[0,1]]\"]", 
                "Nodes Shape: [5, 3]",
                "Edges Count: 14"
                ]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate medium size integer simplex graph (GIS 7 12) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c GIS 7 12")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"GIS\", \"7\", \"12\"]", 
                "Nodes Shape: [18564, 7]",
                "Edges Count: 519792"
                ]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate small integer simplex graph (GIP 3 3), print values to stdout, and check them for corectness":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c GIP 3 3")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"GIP\", \"3\", \"3\"]", 
                "Nodes:",
                "Tensor[system.int] of shape \"[10, 3]\" on backend \"Cpu\"", 
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
                "Neighbors:",
                "@[@[1, 4], @[0, 2, 4, 5], @[1, 3, 5, 6], @[2, 6], @[1, 0, 5, 7], @[2, 1, 4, 6, 7, 8], @[3, 2, 5, 8], @[5, 4, 8, 9], @[6, 5, 7, 9], @[8, 7]]"
                ]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate small fractional simplex graph (GFS 3 3) and print shape to stdout":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c GFS 3 3")
            outputLines = output.splitLines
            reference = @[
                "Running with configuration: @[\"GFS\", \"3\", \"3\"]", 
                "Nodes Shape: [10, 3]",
                "Edges Count: 36"
                ]
        check exitCode == 0
        for i in 0..<reference.len:
            check outputLines[i].stripAnsiEscapeCodes() == reference[i]

    test "generate small fractional simplex graph (GFP 3 3), print values to stdout, and check them for node corectness":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c GFP 3 3")
            outputLines = output.splitLines
            referenceNodes = @[
                @[0.0, 0.0, 1.0],
                @[0.0, 0.333333, 0.666667],
                @[0.0, 0.666667, 0.333333],
                @[0.0, 1.0, 0.0],
                @[0.333333, 0.0, 0.666667],
                @[0.333333, 0.333333, 0.333333],
                @[0.333333, 0.666667, 0.0],
                @[0.666667, 0.0, 0.333333],
                @[0.666667, 0.333333, 0.0],
                @[1.0, 0.0, 0.0]]
    
        check exitCode == 0

        var outputGrid = newSeq[seq[float]](outputLines.len-6)
        for i in 0..<referenceNodes.len:
            for v in outputLines[i+3].replace("|","").splitWhitespace:
                outputGrid[i].add(parseFloat(v))
        
        for i in 0..<referenceNodes.len:
            for j in 0..<3:
                check abs(outputGrid[i][j] - referenceNodes[i][j]) < 0.0001


    test "generate small fractional simplex graph (GFP 3 3), print values to stdout, and check them for neighbor corectness":
        let 
            (output, exitCode) = execCmdEx("./nimplex -c GFP 3 3")
            outputLines = output.splitLines
            referenceNeighbors = @[
                @[1, 4], 
                @[0, 2, 4, 5], 
                @[1, 3, 5, 6], 
                @[2, 6], 
                @[1, 0, 5, 7], 
                @[2, 1, 4, 6, 7, 8],
                @[3, 2, 5, 8], 
                @[5, 4, 8, 9], 
                @[6, 5, 7, 9], 
                @[8, 7]]
        check exitCode == 0
        var outputNeighbors = newSeq[seq[int]](outputLines.len-6)
        for i in 0..<referenceNeighbors.len:
            let parseList = outputLines[14].split("@")[i+2].replace("[", "").replace("],", "").replace("]", "").replace(" ","").split(",")
            for v in parseList:
                outputNeighbors[i].add(parseInt(v))
        
        for i in 0..<referenceNeighbors.len:
            check outputNeighbors[i].toHashSet() == referenceNeighbors[i].toHashSet()

suite "Test NumPy exports corectness for grids":

    test "check if compiled nimplex is present in the current working directory":
        # For Unix systems, the executable is named nimplex, but for Windows it is nimplex.exe
        if hostOS == "windows":
            require fileExists("nimplex.exe")
        else:
            require fileExists("nimplex")

    test "generate auto-named a medium fractional internal simplex grid (IFP 7 11) and export it to NumPy (nimplex_IF_7_11.npy)":
        let (_, exitCode) = execCmdEx("./nimplex -c IFN 7 11")
        check exitCode == 0
        check fileExists("nimplex_IF_7_11.npy")

    test "generate a medium fractional internal simplex grid (IFP 7 11) named testExport.npy and export it to NumPy":
        let (_, exitCode) = execCmdEx("./nimplex -c IFN 7 11 testExport.npy")
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

suite "Test NumPy exports corectness for graphs":

    test "check if compiled nimplex is present in the current working directory":
        # For Unix systems, the executable is named nimplex, but for Windows it is nimplex.exe
        if hostOS == "windows":
            require fileExists("nimplex.exe")
        else:
            require fileExists("nimplex")

    test "generate auto-named a medium fractional simplex graph (GFP 7 11) and export it to NumPy (nimplex_GF_7_11_nodes.npy and nimplex_GF_7_11_neighbors.npy)":
        let (_, exitCode) = execCmdEx("./nimplex -c GFN 7 11")
        check exitCode == 0
        check fileExists("nimplex_GF_7_11_nodes.npy")
        check fileExists("nimplex_GF_7_11_neighbors.npy")

    test "generate a medium fractional simplex graph (GFP 7 11) named testExport.npy and testExport.npy and export it to NumPy":
        let (_, exitCode) = execCmdEx("./nimplex -c GFN 7 11 testExport.npy")
        check exitCode == 0
        check fileExists("testExport_nodes.npy")
        check fileExists("testExport_neighbors.npy")

    test "Verify that hashes of the two NumPy exports are exactly the same binary for nodes":
        let
            hash1 = readFile("nimplex_GF_7_11_nodes.npy").hash.toHex
            hash2 = readFile("testExport_nodes.npy").hash.toHex

        check hash1 == hash2

    test "Verify that hashes of the two NumPy exports are exactly the same binary for neighbors":
        let
            hash1 = readFile("nimplex_GF_7_11_neighbors.npy").hash.toHex
            hash2 = readFile("testExport_neighbors.npy").hash.toHex

        check hash1 == hash2

    test "Verify that the exported NumPy files match the direct result from nimplex (tested for corectness in graph.nim)":
        let 
            result = nimplex.simplex_graph_fractional(7, 11)
            resultNodes = result[0]
            resultNeighbors = result[1]
            loadedNumpyNodes = read_npy[float]("testExport_nodes.npy")
            loadedNumpyNeighbors = read_npy[int]("testExport_neighbors.npy")

        check resultNodes.shape == loadedNumpyNodes.shape
        check resultNodes == loadedNumpyNodes

        var loadedNumpyNeighborsParsed: seq[seq[int]] 
        for nn in loadedNumpyNeighbors.toSeq2D:
            loadedNumpyNeighborsParsed.add(nn.filterIt(it != -1))

        check resultNeighbors.len == loadedNumpyNeighborsParsed.len
        for i in 0..<resultNeighbors.len:
            check resultNeighbors[i].toHashSet() == loadedNumpyNeighborsParsed[i].toHashSet()
    
let t1 = cpuTime()
echo "\n***** CLI BENCHMARK RESULTS *****\n"
echo "Large Graphs:\n" & $initDuration(milliseconds = ((t1 - t0)*1e3).int) & "\n"
echo "-----------------------------------\n"
