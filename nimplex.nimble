# Package

version       = "0.8.0"
author        = "Adam M. Krajewski"
description   = "NIM simPLEX: A concise scientific Nim library (with CLI and Python binding) providing samplings, uniform grids, traversal graphs, and their complexes in compositional (simplex) spaces."
license       = "MIT"
skipDirs      = @["tests", "assets", "docs", "examples"]

# Dependencies
requires "nim >= 2.0.0"
requires "arraymancer >= 0.7.33"
requires "nimpy"