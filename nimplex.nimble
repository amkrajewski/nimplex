# Package

version       = "0.5.0"
author        = "Adam M. Krajewski"
description   = "NIM simPLEX: A concise scientific Nim library (with CLI and Python binding) providing samplings, uniform grids, and traversal graphs in compositional (simplex) spaces."
license       = "MIT"
srcDir        = "src"
skipDirs      = @["tests", "assets", "docs", "examples"]

# Dependencies
requires "nim >= 2.0.0"
requires "arraymancer >= 0.7.3"
requires "nimpy"
