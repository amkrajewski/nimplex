## An example of using Python library (pycalphad) from Nim to calculate equilibrium phases of Ti40 W50 Zr10.
## It is meant to be called from the `examples` directory. No special flags are needed to run this script,
## as long as the `myPycalphadCallable.py` file is in the working directory, dependencies are installed, and
## the TDB files are in the correct location. To run, simply `nim c -r -d:release extras/pycalphadFromNim.nim`.

import nimpy
import std/os
import std/times

let sys = pyImport("sys")
discard sys.path.append(getCurrentDir())

let pycalphadCap = pyimport("myPycalphadCallable")

let t0 = cpuTime()

# elementalSpaceComponents = ["Ti", "Zr", "Hf", "W", "Nb", "Ta", "Mo"]
let eqPhases = pycalphadCap.equilibrium_callable([0.4, 0.1, 0, 0.5, 0, 0, 0])

echo "Equilibrium phases of Ti40 W50 Zr10:"
for p in eqPhases:
  echo p

let t1 = cpuTime()
echo "\n", "First Calculation Time: ", $initDuration(milliseconds = ((t1 - t0)*1e3).int)

for _ in 0..10:
  discard pycalphadCap.equilibrium_callable([0.5, 0, 0, 0.5, 0, 0, 0])
let t2 = cpuTime()
echo "Serial Calculation Average Time: ", $initDuration(milliseconds = ((t2 - t1)*1e3/10).int), "\n"

