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

func findSubspace[T](
        lst: seq[T],
        maxDim: int = 3
    ): seq[seq[T]] =
    ## Finds all possible subspaces of the input list `lst` up to the dimension `maxDim`. The list can be of any type and macro
    ## will handle it returning a list of lists of the same type. The default `maxDim` is set to 3 because for larger spaces, especially
    ## high dimensional ones, the number of subspaces grows rapidly (e.g. 1,956 for d=6) and computer memory can be exhausted quicker than
    ## anticipated when working on it in more intense steps; user can set it to any value they want though.
    let n = lst.len
    var sublist = newSeq[T]()
    
    for i in 1 ..< (1 shl n):
        sublist = @[]
        for j in 0 ..< n:
            if (i and (1 shl j)) != 0:
                sublist.add(lst[j])
        if sublist.len <= maxDim:
            result.add(sublist)


func isSubspace(subSpace: seq[int], space: seq[int]): bool =
    ## Small helper function to check if the input `subSpace` is a subspace of the `space`.
    for i in subspace:
        if i notin space:
            return false
    return true


func compareNodes(
        a: int, 
        b: int, 
        nodeCoordinateTensor: Tensor[int], 
        priorityList: seq[int]
    ): int =
    ## Comparator function to evaluate the order of two node numbers `a` and `b` based on their coordinates in the compositional space given in 
    ## `nodeCoordinateTensor` of `int`s (e.g. `simplex_grid` output or the first element in the `simplex_graph` output tuple), ranked by
    ## the `priorityList` establishing order of dimensions in the ordered subspace one may be interested in. 
    let 
        p1 = nodeCoordinateTensor[a, _].squeeze()
        p2 = nodeCoordinateTensor[b, _].squeeze()

    for p in priorityList:
        if p1[p] > p2[p]:
            return -1
        if p1[p] < p2[p]:
            return 1
    return 0


func sortNodes(
        nodes: var seq[int],
        coordinateTensor: Tensor[int], 
        priorityList: seq[int]
    ): void =
    ## Sorts a `var` `seq` of node numbers `nodes` using the `compareNodes` function.
    nodes.sort((x, y) => compareNodes(x, y, coordinateTensor, priorityList))


func permutations[T](s: seq[T]): seq[seq[T]] =
    ## Generates all possible permutations of the input sequence `s` using a recursive algorithm. Can accept any type of sequence thanks
    ## to Nim's powerful type system and reurns a sequence of sequences of the same type.
    if s.len == 0:
        return @[]
    if s.len == 1:
        return @[s]
    
    var xs: seq[T]
    for i in 0..<s.len:
        let x = s[i]
        xs = s[0..<i] & s[i+1..^1]
        for p in permutations(xs):
            result.add(@[x] & p)