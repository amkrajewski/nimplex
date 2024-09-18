import std/tables
import std/math
import std/algorithm
import std/sugar
import arraymancer/Tensor
import ../nimplex

## This submodule contains utility functions related to stitching of the compositional graphs to form graph complexes, so that
## much more complex graphs can be built from simpler ones while retaining homogeneous structure of the space. Furthermore,
## one can keep track of provenance of the subgraphs and use this information to deploy computational (e.g., ML) models
## on per-subgraph basis, which should be extremely useful for (a) combining the power of many specialized models and (b)
## creating stacked spaces for multi-step problems broken down into individual steps.
## 
## **Navigation:** [nimplex](../nimplex.html) (core library) | [docs/changelog](../docs/changelog.html) | [utils/plotting](plotting.html)
## 

func generateAlphabetSequence(length: int): seq[string] =
    assert length <= 52, "The number of characters requested to name system components exceeds the number of available symbols (52) - uppercase and lowercase letters)."
    const alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    for i in 0..<length:
        result.add($alphabet[i])

