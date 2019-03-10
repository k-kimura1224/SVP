# Branch-and-Bound Solver for Shortest Vector Problems

## Compile
1. Write `CXX` (for C++14), `CXXFLAGS` and `CXXLDFLAGS` (for BLAS) in Makefile.
2. Execute `make OPT=opt` to compile this solver.

## Usage
`SVPSOLVER.opt -f filename [-p nthreads] [-t timelimit] [-m memory] [-H AF] [-q] [-h]`

- -f filename: filename of dat-file

## Option
- -p nthreads: number of threads (default: 1)
- -t timelimit: timelimit for solving (default: 5000(s))
- -m memory: max value for memory (default: 8(GB))
- -H AF: approximation factor to execute heuristic methods (default: 0.95)  
a
- -q: quiet mode (default: false)
- -h: display usage

## Example
- `./bin/SVPSOLVER.opt -f sample.dat`
- `./bin/SVPSOLVER.opt -f sample.dat -p 4 -m 8 -t 5000`
- `./bin/SVPSOLVER.opt -f sample.dat -p 32 -m 128 -t 86400 -H 1.04`
