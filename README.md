# nimplex

NIM simPLEX: A concise scientific library (with CLI) providing uniform density grid and sampling in simplex spaces.

[![(Linux) Grid, Graph, and CLI Tests](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Linux.yaml/badge.svg)](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Linux.yaml)

[![(MaxOS) Grid, Graph, and CLI Tests](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Apple.yaml/badge.svg)](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Apple.yaml)

[![(Windows) Grid, Graph, and CLI Tests](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Windows.yaml/badge.svg)](https://github.com/amkrajewski/nimplex/actions/workflows/testingOnPush_Windows.yaml)


## Installation
There are several **easy** ways to quickly get *nimplex* up and running on your system. The choice depends primarily on your preffered way of interacting with the library (CLI, Nim, or Python) and your system configuration.

If you happen to be on one of the common systems (for which we auto-compile the binaries) and you do not need to modify anything in the source code, you can simply download the latest release from the [nimplex GitHub repository](https://github.com/amkrajewski/nimplex)
and run the executable (*nimplex* / *nimplex.exe*) or Python library (*nimplex.so* / *nimplex.pyd*) directly just by placing it in your working directory and using it as:

1. An **interactive command line interface (CLI) tool**, which will guide you through how to use it if you run it without any arguments like so (on Linux/MacOS):
   ```cmd
   ./nimplex   
   ```
   or with a concise configuration defining the task type and parameters (explained later in [Usage in Nim](#usage-in-nim)):
   ```cmd
   ./nimplex -c IFP 3 10
   ```
2. An **compiled Python library**, which you can import and use in your Python code like so:
   ```python
   import nimplex
   ```
   and immediately use the functions provided by the library, as described in [Usage in Python](#usage-in-python):
   ```python
   nimplex.simplex_internal_grid_fractional(dim=3, ndiv=10)
   ```

If the above doesn't work for you, or you want to modify the source code, you can compile the library yourself fairly easily in a couple minutes. 
The only requirement is to have [Nim](https://nim-lang.org/) installed on your system
([Installation Instructions](https://nim-lang.org/install.html)) which can be done on most Linux distributions with a single command:
```cmd
apt-get install nim
```
or on MacOS, assuming you have [Homebrew](https://brew.sh/) installed:
```cmd
brew install nim
```

Then, you can use the boundeled [Nimble](https://github.com/nim-lang/nimble) tool (pip-like package manager for Nim) to install two top-level dependencies: 
[arraymancer](https://github.com/mratsim/Arraymancer), which is a powerful N-dimensional array library, and [nimpy](https://github.com/yglukhov/nimpy) which 
helps with the Python bindings. You can do it with a single command:
```cmd
nimble install  -y arraymancer nimpy
```

Finally, you can clone the repository and compile the library with:
```cmd
git clone https://github.com/amkrajewski/nimplex
cd nimplex
nim c -r -d:release nimplex.nim -benchmark
```
which will compile the library and run a few benchmarks to make sure everything runs smoothly. You should then see a compiled binary file `nimplex` in the current directory which exposes the CLI tool.
If you want to use the Python bindings, you can compile the library with slightly different flags (depending on your system configuration) like so for Linux/MacOS:
```cmd
nim c --d:release --threads:on --app:lib --out:nimplex.so nimplex
```
and you should see a compiled library file `nimplex.so` in the current directory which can be immediately imported and used in Python.