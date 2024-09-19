import std/tables
import std/math
import std/algorithm
import std/sugar
import std/times
import std/strformat
import arraymancer/Tensor
import ../nimplex

## This submodule contains utility functions related to stitching of the compositional graphs to form graph complexes, so that
## much more complex graphs can be built from simpler ones while retaining homogeneous structure of the space. Furthermore,
## one can keep track of provenance of the subgraphs and use this information to deploy computational (e.g., ML) models
## on per-subgraph basis, which should be extremely useful for (a) combining the power of many specialized models and (b)
## creating stacked spaces for multi-step problems broken down into individual steps.
## 
## **Navigation:** [nimplex](../nimplex.html) (core library) | [docs/changelog](../docs/changelog.html) | [utils/plotting](plotting.html) | [utils/stitching](stitching.html)
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
    ## Finds all possible (unordered) subspaces of the input list `lst` up to the dimension `maxDim`. The list can be of any type and macro
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


proc findStitchingPoints*(
    dim: int, 
    ndiv: int,
    maxDim: int = 3,
    components: seq[string] = generateAlphabetSequence(dim)
        ): Table[string, seq[int]] =
    let L: int = binom(ndiv+dim-1, dim-1)
    var 
        x = zeros[int](dim)
        stitchTableInt = initTable[seq[int], seq[int]]()
        maxSys: seq[int] = @[]
        grid: Tensor[int] = newTensor[int]([L, dim])

    # Generate a space vector with all components present
    for i in 0..<dim:
        maxSys.add(i)

    # Generate all (un-ordered) subspaces of the space up to the given dimension
    for subSys in findSubspace(maxSys, maxDim):
        stitchTableInt[subSys] = @[]

    # Start the generation of compositional grid that will be used for sorting purposes
    x[dim-1] = ndiv
    for j in 0..<dim:
        grid[0, j] = x[j]

    # Start the generation of the stitch table
    for sys in stitchTableInt.keys:
        if nonZeroComps(x).isSubspace(sys):
            stitchTableInt[sys].add(0)

    # Iterate over all possible compositions of the space using NEXCOM algorithm (see the nimplex paper for details)
    # while building the stitching table through assigning current point to the appropriate subspace and its subspaces
    var h = dim
    for i in 1..<L:
        h -= 1
        let val = x[h]
        x[h] = 0
        x[dim-1] = val - 1
        x[h-1] += 1
        for j in 0..<dim:
            grid[i, j] = x[j]
        if val != 1:
            h = dim
        for sys in stitchTableInt.keys:
            if nonZeroComps(x).isSubspace(sys):
                stitchTableInt[sys].add(i)
    
    # For every unordered subspace we found, find all of its permutations and sort the associated nodes based on the 
    # order of elements in said permutations
    for sys in stitchTableInt.keys:
        var sortedSys: seq[int] = stitchTableInt[sys]
        let permutations = sys.permutations
        for p in permutations:
            sortedSys.sortNodes(grid, p)
            result[space2name(p, components)] = sortedSys

if isMainModule:
    let stitch = findStitchingPoints(5, 5, 3)
    for space in stitch.keys:
        echo fmt"{space:<7} -> {stitch[space]}"
    
    echo "\n\n"

    let stitch2 = findStitchingPoints(3, 4, 3, @["Ti", "V", "Cr"])
    for space in stitch2.keys:
        echo fmt"{space:<7} -> {stitch2[space]}"

    let t0 = cpuTime()
    let stitch3 = findStitchingPoints(6, 9, 4)
    let t1 = cpuTime()
    var 
        stitchCount: int = 0
        subspaceCount: int = 0
    for space in stitch3.keys:
        subspaceCount += 1
        stitchCount += stitch3[space].len
    echo "\n", fmt"Benchmark: Found {stitchCount} stitch points on all permutations of {subspaceCount} quaternary, ternary, binary, and unary subspaces of a 6-dimensional space with 9 divisions per dimension in {t1-t0} seconds."