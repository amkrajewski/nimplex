import arraymancer/Tensor

## This submodule contains utility functions related to plotting of the compositional data.

proc simplex2cartesian*(simplexPoints: Tensor[float]): Tensor[float] =
    ## Converts Arraymancer `Tensor[float]` of `simplexPoints` with fractional coordinates (e.g., from grid or random sampling) to points in a cartesian space (within unit n-sphere)
    ## for purposes of plotting in the much more common cartesian space.
    let dim = simplexPoints.shape[1]
    assert dim > 0, "0-component simplex is not defined"
    
    case dim:
        of 1:
            # 1-component simplex is a point, so all are at the origin (0, 0, 0, ...)
            return zeros[float](simplexPoints.shape)
        of 2:
            # 2-component simplex is a line segment of lenght 1, so points are just either the first or second (mirrored) column of the input Tensor
            return simplexPoints[_, 0]
        of 3:
            # 3-component simplex can be projected onto the 2D plane as an equilateral triangle inside a unit circle, so we can use the following formula to convert from barycentric to cartesian coordinates:
            let transformationTensor = @[
                @[0.0, 1.0],
                @[0.866025, -0.5],
                @[-0.866025, -0.5]
                ].toTensor()
            return simplexPoints * transformationTensor
        of 4:
            # 4-component simplex can be projected onto the 3D space as a tetrahedron inside a unit sphere, so we can use the following formula to convert from barycentric to cartesian coordinates:
            # {{0.942809, 0., -0.333333}, {-0.471405, 0.816497, -0.333333}, {-0.471405, -0.816497, -0.333333}, {0., 0., 1.}}
            let transformationTensor = @[
                @[0.94280904, 0.0, -0.3333334],
                @[-0.47140452, 0.8164966, -0.3333334],
                @[-0.47140452, -0.8164966, -0.3333334],
                @[0.0, 0.0, 1.0]
                ].toTensor()
            return simplexPoints * transformationTensor
        
        else:
            raise newException(ValueError, "(input dim" & $dim & ") The simplex2cartesian utility function is defined only for 1, 2, 3, and 4-component simplexes corresponding to 0D, 1D, 2D, and 3D spaces respectively, as the main motivation for this function is to visualize the simplex in 2D or 3D space. If you need to convert a simplex of higher dimensionality, please open an issue on GitHub or implement it yourself and submit a PR.")

# PYTHON BINDINGS
when not defined(nimdoc):
    import nimpy
    
    proc cartesian2simplex_py*(simplexPoints: seq[seq[float]]): seq[seq[float]] {.exportpy.} =
        return simplex2cartesian(simplexPoints.toTensor()).toSeq2D()