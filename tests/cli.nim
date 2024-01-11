import std/[os, osproc]
import std/unittest
import std/strutils

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
