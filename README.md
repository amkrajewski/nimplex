# nimplex

[![Static Badge](https://img.shields.io/badge/API-Documentation-e?logo=bookstack&logoColor=white&link=https%3A%2F%2Famkrajewski.github.io%2Fnimplex%2F)](https://nimplex.phaseslab.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<a href="https://nimble.directory/pkg/nimplex"><img src="https://raw.githubusercontent.com/amkrajewski/nimble-badge/master/badges/nimplex.svg" alt="Nimble Package" height="20"></a>
[![Arxiv](https://img.shields.io/badge/arXiv-2402.03528-8F1515?style=flat&logo=arxiv&logoColor=red)](https://doi.org/10.48550/arXiv.2402.03528) 
[![npj Unconventional Computing](https://img.shields.io/badge/npj_Unconventional_Computing-10.1038%2Fs44335--024--00012--2-blue)](https://doi.org/10.1038/s44335-024-00012-2)

[![(Linux) Grid, Graph, and CLI Tests](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Linux.yaml/badge.svg)](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Linux.yaml)
[![(MaxOS) Grid, Graph, and CLI Tests](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Apple.yaml/badge.svg)](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Apple.yaml)
[![(Windows) Grid, Graph, and CLI Tests](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Windows.yaml/badge.svg)](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Windows.yaml)

**NIM** sim**PLEX**: A concise high-performance scientific Nim library (with CLI and Python binding) providing samplings, uniform grids, traversal graphs, and more in compositional (simplex) spaces, where traditional methods designed for euclidean spaces fail or otherwise become impractical.

Such spaces are considered when an entity can be split into a set of distinct components (a composition), and they play a critical role in many disciplines of science, engineering, and mathematics. For instance, in materials science, chemical composition refers to the way a material (or, more generally, matter) is split into distinct components, such as chemical elements, based on considerations such as fraction of atoms, occupied volume, or contributed mass. And in economics, portfolio composition may refer to how finite capital is split across assets, such as cash, equity instruments, real estate, and commodities, based on their monetary value.

## Quick Start

If you have a GitHub account, you can get started with `nimplex` very quickly by just clicking the button below to launch a CodeSpaces environment with everything installed (per instructions in [Reproducible Installation](#reproducible-installation) section) and ready to go! From there, you can either use the CLI tool (as explained in [CLI](#cli) section) or import the library in Python (as explained in [Usage in Python](#usage-in-python) section) and start using it right away. Of course, it also comes with a full Nim compiler and VSCode IDE extensions for Nim, so you can efortlessely modify/extend the source code and re-compile it if you wish.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/amkrajewski/nimplex?quickstart=1)

## Installation
There are several **easy** ways to quickly get *nimplex* up and running on your system. The choice depends primarily on your preffered way of interacting with the library (CLI, Nim, or Python) and your system configuration.

### Reproducible Installation (recommended)

The recommended way is compiling the library yourself, which may sound scary but is fairly easily and the whole process should not take more than a couple of minutes. 

#### Nim (compiler)

First, you need to install [Nim](https://nim-lang.org/) language compiler which on most **Unix** (Linux/MacOS) systems is very straightforward.

- On **MacOS**, assuming you have [Homebrew](https://brew.sh/) installed, simply:
  ```sh
  brew install nim
  ```

- Using **`conda`**, [`miniconda`](https://docs.anaconda.com/miniconda/), `mamba`, or [`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) cross-platform package manager:
  ```sh
  conda install -c conda-forge nim
  ```

- On most **Linux** distributions, you should also be able to use your built-in package manager like `pacman`, `apt`, `yum`, or `rpm`; however, the default channel/repository, especially on enterprise systems, may have an unsupported version (`nim<2.0`). While we [do test `nimCSO` with `1.6` versions too](https://github.com/amkrajewski/nimCSO/blob/main/.github/workflows/testingOnPush_Linux.yaml), your experience may be degraded, so you may want to update it or go with another option.

- You can, of course, also build it yourself from [`nim` source code](https://github.com/nim-lang/Nim)! It is relatively straightforward and fast compared to many other languages. 

On **Windows**, you may consider using [`WSL`](https://learn.microsoft.com/en-us/windows/wsl/about), i.e., Windows Subsystem for Linux, which is strongly recommended, interplays well with VS Code, and will let you act as if you were on Linux. If you need to use Windows directly, you can follow these [installation instructions](https://nim-lang.org/install_windows.html).

#### `nimplex`

Then, you can use the bundled [Nimble](https://github.com/nim-lang/nimble) tool (package manager for Nim, similar to Rust's `crate` or Python's `pip`) to install two top-level `nim` dependencies:
- [`arraymancer`](https://github.com/mratsim/Arraymancer), which is a powerful N-dimensional array library
- [`nimpy`](https://github.com/yglukhov/nimpy) which helps with the Python bindings. 

It's a single command:

```sh
nimble install --depsOnly
```

or, explicitly:

```sh
nimble install -y arraymancer nimpy
```

Finally, you can clone the repository and compile `nimplex` with:
```cmd
git clone https://github.com/amkrajewski/nimplex
cd nimplex
nim c -r -d:release nimplex.nim --benchmark
```
which will compile the library and run a few benchmarks to make sure everything runs smoothly. You should then see a compiled binary file `nimplex` in the current directory which exposes the CLI tool.

If you want to use the **Python bindings**, you can now compile the library with slightly different flags (*depending on your system configuration*) like so for Linux/MacOS:
```cmd
nim c --d:release --threads:on --app:lib --out:nimplex.so nimplex
```
and you should see a compiled library file `nimplex.so` in the current directory which can be immediately imported and used in Python as explained later. For Windows and other platforms, consult [`nimpy`](https://github.com/yglukhov/nimpy) documentation on what flags and formats should be used.

### Pre-Compiled Binaries (quick but not recommended)

If you happen to be on one of the common systems (for which we auto-compile the binaries) and you do not need to modify anything in the source code, there is a good chance you can simply download the latest release from the [nimplex GitHub repository](https://github.com/amkrajewski/nimplex)
and run the executable (*nimplex* / *nimplex.exe*) or Python library (*nimplex.so* / *nimplex.pyd*) directly just by placing it in your working directory and using it as:

1. An **interactive command line interface (CLI) tool**, which will guide you through how to use it if you run it without any arguments like so (on Linux/MacOS):
   ```cmd
   ./nimplex   
   ```
   or with a concise configuration defining the task type and parameters (explained later in [Usage in Nim](#usage-in-nim)):
   ```cmd
   ./nimplex -c IFP 3 10
   ```
2. An **compiled Python library for Unix**, which you can import and use in your Python code like so:
   ```python
   import nimplex
   ```
   and immediately use the functions provided by the library, as described in [Usage in Python](#usage-in-python):
   ```python
   nimplex.simplex_internal_grid_fractional(dim=3, ndiv=10)
   ```

## Capabilities
***Note:*** Full technical discussion of methods and motivations is provided in the manuscript. The sections below are meant to provide a concise overview of the library's capabilities.

The library provides a growing number of methods specific to compositional (simplex) spaces:
1. **Monte Carlo sampling** is the simplest method conceptually, where points are rendomly sampled from a simplex. In low dimensional cases, this can be accomplished by sampling from a uniform distribution in (d-1)-Cartesian space and then rejecting points outside the simplex (left panel below). However, in this approach, the inefficiency growth is **factorial** with the dimensionality of the simplex space. Instead, some try to sample from a uniform distribution in (d)-Cartesian space and normalize the points to sum to 1, however, this leads to over-sampling in the center of each simplex dimension (middle panel below). 

   One can, however, fairly easily sample from a special case of Dirichlet distribution, as explained in the manuscript, which leads to uniform sampling in the simplex space (right panel below). **Nimplex can sample around 10M points per second in 9-dimensional space** on a modern CPU.

   <img src="https://raw.githubusercontent.com/amkrajewski/nimplex/main/assets/Fig1.png" alt="Monte Carlo Sampling" width="800" style="display: block; margin-left: auto; margin-right: auto;"/>
   
2. **Simplex / Compositional Grids** are a more structured approach to sampling, where all possible compositions quantized to a given resolution, like 1% for 100 divisions per dimension, are generated. This is useful for example when one wants to map a function over the simplex space. In total `N_S(d, n_d) = \binom{d-1+n_d}{d-1} = \binom{d-1+n_d}{n_d}` are generated, where `d` is the dimensionality of the simplex space and `n_d` is the number of divisions per dimension. Nimplex uses a modified version of NEXCOM algorithm to do that procedurally (see manuscript for details) and can generate around **5M points per second in 9-dimensional space** on a modern CPU. A choice is given between generating the gird as a list of **integer** numbers of quantum units (left panel below) or as a list of **fractional positions** (right panel below). 
   <img src="https://raw.githubusercontent.com/amkrajewski/nimplex/main/assets/Fig2.png" alt="Integer and Fractional Simplex Grids in Ternary Space" width="800" style="display: block; margin-left: auto; margin-right: auto;"/>

3. **Internal Simplex / Compositional Grids** are a modification of the above method, where only points inside the simplex, i.e. all components are present, are generated. This is useful in cases where, one cannot discard any component entirely, for instance, because manufacturing setup has minimum feed rate (leakage). Nimplex introduces a new algorithm to generate these points procedurally (see manuscript for details) based on further modification of NEXCOM algorithm. 

   In total `N_I(d, n_d) = \binom{n_d-1}{d-1}` are generated, critically without any performance penalty compared to the full grid, which can reach orders of magnitude when `d` approaches `n_d`. Similar to the full grid, a choice is given between generating the gird as a list of **integer** numbers of quantum units or as a list of **fractional positions**.

4. **Simplex / Compositional Graphs** generation is ***the most critical capability***, first introduced in the nimplex manuscript. They are created by using combinatorics and disocvered patterns to assign edges between all neighboring nodes during the simplex grid (graph nodes) generation process. Effectively, a traversal graph is generated, spanning all possible compositions (given a resolution) creating an extremely efficient representation of the problem space, which allows deployment of numerous graph algorithms. 

   <img src="https://raw.githubusercontent.com/amkrajewski/nimplex/main/assets/Fig3.png" alt="Simplex Graph for Ternary Space" width="800" style="display: block; margin-left: auto; margin-right: auto;"/>
   
   Critically, unlike the O(N^2) distance-based graph generation methods, this approach **scales linearly** with the resulting number of nodes. Because of that, it is extremely efficient even in high-dimensional spaces, where the number of edges goes into trillions and beyond. Nimplex can **both generate and find neighbors** for around **2M points per second in 9-dimensional space** on a modern CPU. 

5. **Limited Simplex / Compositional Graphs** where each component (i.e. simplex dimension) can be individually limited by minimum and maximum values, as depicted in the figure below. This allows for efortless creation of graphs representing asymmetric subspaces corresponding to, for instance, “all HEA with at least 5% and at most 55% of each of 8 metallic components but less than 5% of boron and between 4% and 12% of carbon”.

   <img src="https://raw.githubusercontent.com/amkrajewski/nimplex/main/assets/small_GL.png" alt="Limited Simplex Graph for Ternary Space" width="800" style="display: block; margin-left: auto; margin-right: auto;"/>
   
6. **Graph Complexes**, combining multiple individual compositional graphs can also be constructed through the `utils/stitching` module which computes fixed-orientation subspaces that can be joined together,
   e.g., A-B-C in A-B-C-D with A-B-C in F-A-C-E-B. As explored in [our npj Unconventional Computing article](https://doi.org/10.1038/s44335-024-00012-2), such combined representations, allowing for even different dimensionalities of joined spaces, can can then be used to efficeintly encode complex problem spaces where some prior assumptions and knowledge are available. In the Example #2 from our article, inspired by problem of joining titanium with stainless steel in [10.1016/j.addma.2022.102649](https://doi.org/10.1016/j.addma.2022.102649), using 3-component spaces, one can encode 3 separate paths where some components are shared in predetermined fashion. This efficiently encodes the problem space in form of a single simplex graph complex (right panel below) with a single consistent structure, that can be directly used for pathfinding and other graph algorithms just like any other graph.

   <img src="https://raw.githubusercontent.com/amkrajewski/nimplex/main/assets/Fig4.png" alt="Simplex Graph Complex" width="1000" style="display: block; margin-left: auto; margin-right: auto;"/>

   With such graph representation, one can very easily deploy any scientific library for graph exploration, constrained and biased by models operating in the elemental space mapping `nimplex` provides. A neat and concise demonstration of this is provided in the [`02.AdditiveManufacturingPathPlanning.ipynb`](examples/02.AdditiveManufacturingPathPlanning.ipynb) under `examples` directory, where thermodynamic phase stability models constrain a 4-component (tetrahedral) design space existing in 7-component chemical space and property model related to yield strength (RMSAD) is used to bias designed paths towards objectives like property maximization or gradient minimization with extremely concise code simply modifying the weights on unidirectional edges in the graph. For instance, the figure below (approximately) depicts the shortest path through a subset of tetrahedron formed by solid solution phases, later stretched in space proportionally to RMDAS gradient magnitude.

   <img src="https://raw.githubusercontent.com/amkrajewski/nimplex/main/assets/Fig5.png" alt="Gradient Magnitude Stretched Graph with Shortest Path" width="1000" style="display: block; margin-left: auto; margin-right: auto;"/>

   A real-world example of using such context can be found in [our paper on AMMap tool](https://doi.org/10.31224/4453), where a graph complex consisting of all possible ternary (3-component) simplex graphs in a 7-component space are joined together, representing a continous space of all possible ternary compositions of 7 elements. In the paper, it is then used to deploy thermodynamic calculations through a search algorithm, with aim of finding the shortest path between two compositions of interest. The resulting graph complex with calculation results overlaid is shown below.

   <img src="https://raw.githubusercontent.com/amkrajewski/nimplex/main/assets/Fig6.png" alt="AMMap Simplex Graph Complex" width="500" style="display: block; margin-left: auto; margin-right: auto;"/>

> Several other methods are in testing and will likely be added in the future releases. If you have any suggestions, please open an issue on GitHub as we are always soliciting new ideas and use cases based on real-world problems in the 
scientific computing community.

## Usage in Nim
Usage within Nim is fairly straightforward. You can install it using Nimble as explained earlier, or install it directly from GitHub, making sure to use the slightly modified `@#nimble` branch:
```cmd
nimble install -y https://github.com/amkrajewski/nimplex@#nimble
```
or, if you wish to modify the source code, you can simply download the core file `nimplex.nim` and place it in your own code, as long as you have the dependencies installed, since it is standalone. 
**Then simply follow the [API documentation (amkrajewski.github.io/nimplex)](https://amkrajewski.github.io/nimplex/) which goes over all core functions and extra utilities like [nimplex/utils/plotting](https://amkrajewski.github.io/nimplex/utils/plotting.html) and [nimplex/utils/stitching](https://amkrajewski.github.io/nimplex/utils/stitching.html).**

## Usage in Python
To use the library in Python, you can interact with it just like any other Python library. All input/output types are native Python types, so no additional conversion is necessary!. Once you have the library installed and imported,
**simply follow the [API documentation](https://amkrajewski.github.io/nimplex/), with an exception that you need to add `_py` to the function names.** If you happen to forget adding `_py`, the Python interpreter will throw an error with a suggestion to do so. A couple of additional conveninece functions are listed under [nimplex/#usage-in-python](https://amkrajewski.github.io/nimplex/#usage-in-python). 

## CLI

### Interactive
Using Nimplex through the CLI relies on the same core library, but provides a simple interface for users who do not want to write any code. It can be used interactively, where the user is guided through the configuration process by 
just running the executable without any arguments:
```cmd
./nimplex
```

### Configured
Or it can be run with a concise configuration defining the task type and parameters. The configuration is a 3-letter string and 2-3 additional parameters, as explained below.
- **3-letter configuration**: 
    1. Grid type or uniform random sampling:
        - **F**: Full grid (including the simplex boundary)
        - **I**: Internal grid (only points inside the simplex)
        - **R**: Random/Monte Carlo uniform sampling over simplex
        - **G**: Graph (list of grid nodes and list of their neighbors)
        - **L**: Limited graph (graph with integer or fractional limits imposed on the compositions)
    2. Fractional or Integer positions:
        - **F**: Fractional grid/graph (points are normalized to fractions of 1)
        - **I**: Integer grid/graph (points are integers)
    3. Print full result, its shape, or persist in a file:
        - **P**: Print (presents full result as a table)
        - **S**: Shape (only the shape / size information)
        - **N**: Persist to NumPy array file ("nimplex_<configFlags>.npy" or optionally a custom path as an additonal argument)
- **Simplex Dimensions / N of Components**: An integer number of components in the simplex space.
- **N Divisions per Dimension / N of Samples**: An integer number of either:
    1. Divisions per each simplex dimension for grid or graph tasks (F/I/G__)
    2. Number of samples for random sampling tasks (R__)
- **(for limited graphs only) Limits**: A string list of pairs of either integers or floats, depending on the grid type specified by the second letter of 
  the configuration string. Each pair represents the minimum and maximum values for each simplex dimension. The number of pairs must match the 
  number of simplex dimensions. You can use `[`, `{`, or `@[` for the list. E.g., `[[0, 1],[0,0.555],[0.222,0.888]]` for fractional grid limits
  or `{{0,12},{2,10},{3,9}}` for integer grid limits. Setting the limits to less than 0 or greater than 1/nDiv will pass with a warning, but will
  not throw an error. The integer limits are exactly inclusive, while the fractional limits inclusive of the *closeset* grid point, i.e. max limit of 
  0.660 will let 1/3=0.667 pass for low resolution grids, but not high resolution ones.
- **(optional) NumPy Array Output Filename**: A custom path to the output NumPy array file (only for __N tasks). 

For instance, to generate a 3-dimensional internal fractional grid with 10 divisions per dimension and persist it to a NumPy array file, you can run:
```cmd
./nimplex -c IFN 3 10
```
and the output will be saved to `nimplex_IF_3_10.npy` in the current directory. If you want to save it to a different path, you can provide it as an additional argument:
```cmd
./nimplex -c IFN 3 10 path/to/outfile.npy
```
Or if you want to print the full result to the console, allowing you to pipe it to virtually any other language or tool as plain text, you can run:
```cmd
./nimplex -c IFP 3 10
```

### Auxiliary Flags
You can also utilize the following auxiliary flags:
- `--help` or `-h` --> Show help.
- `--benchmark` or `-b` --> Run a set of tasks to benchmark performnace (`simplex_grid(9, 12)`, `simplex_internal_grid(9, 12)`, `simplex_sampling_mc(9, 1_000_000)`, `simplex_graph(9, 12)`) and
    compare performance across implementations (`simplex_graph(3, 1000)` vs `simplex_graph_3C(1000)`).
