# Package

version       = "0.7.1"
author        = "Adam M. Krajewski"
description   = "NIM simPLEX: A concise scientific Nim library (with CLI and Python binding) providing samplings, uniform grids, traversal graphs, and their complexes in compositional (simplex) spaces."
license       = "MIT"
srcDir        = "src"
skipDirs      = @["tests", "assets", "docs", "examples"]
bin           = @["nimplex"]
installExt    = @["nim"]

# Dependencies
requires "nim >= 2.0.0"
requires "nimcuda >= 0.1.4 & <= 0.1.9"
requires "arraymancer >= 0.7.3 & <= 0.7.32"
requires "nimpy"
