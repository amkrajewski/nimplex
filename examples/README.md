# nimplex Examples 

## Quick Start

This is a Quick Start guide to get you up and running with the nimplex software. If you are running it in the Codespaces environment, you should now open the 
[**`01.QuickStart.ipynb`**](01.QuickStart.ipynb) notebook inside this directory and follow the instructions there. **Please note, once propted to select a kernel (after running the first code cell), you should select the `base` kernel (which should show up as recommended).**

## Additive Manufacturing Path Planning Made Efortless

**In the [`02.AdditiveManufacturingPathPlanning.ipynb`](02.AdditiveManufacturingPathPlanning.ipynb) tutorial, we will demonstrate how effortless it is to dramatically speed up the exploration of feasible compositional spaces in high dimensional spaces through employing `nimplex`'s graph representations that abstract the underlying problem and dimensionality.**

**We will also design several neat, mathematically optimal (given some criteria) paths in a 7-component chemical space connecting two alloys of interest by mixing 4 fixed-composition alloy powders to create a tetrahedral attainable/design space. The beauty of this approach is that at no point (except for plotting in 3D for "human consumption") will we explicitly consider the dimensionality or the distance as the connectivity between the points in the space has been abstracted into graph adjacency. If you wish to add another alloy to the design process, you add it to the list, and you are done :)**

In addition to leveraging the graph representation to glide around phase boundaries of infeasible regions, we will also demonstrate path planning for:
- Shortest Path
- Property Gradient Minimization Bias Path
- Property Value Maximization Bias Path
- (and several more mentioned as bonus exercises!)
