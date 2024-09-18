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
    ## Generates a sequence of strings of length `length` containing all UPPER and lower case letters of the alphabet. 
    ## The primary purpose of this function is to generate unique names for space (A-B-C-D) and its ordered subspaces
    ## (e.g., A-B, B-A) if no custom names are provided by the user.
    assert length <= 52, "The number of characters requested to name space components exceeds the number of available symbols (52) - uppercase and lowercase letters)."
    const alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    for i in 0..<length:
        result.add($alphabet[i])


func nonZeroComps(p: Tensor[int]): seq[int] =
    ## Finds the indices of non-zero components in the input tensor `p` to determine which (unordered) space/subspace it belongs to.
    for i in 0..<p.shape[0]:
        if p[i] != 0:
            result.add(i)


func space2name(space: seq[int], components: seq[string]): string =
    ## Converts the space/subspace represented by the indices in `space` to a string representation using the names of the `components`.
    ## It can be used in conjunction with `generateAlphabetSequence` too.
    for i in space:
        result &= components[i]
        if i != space[space.len-1]:
            result &= "-"