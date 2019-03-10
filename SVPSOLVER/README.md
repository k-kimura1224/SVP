# Branch-and-Bound Solver for Shortest Vector Problems

## Compile
Write ***`CXX` Makefile path of the directory of your scip in SCIPDIR
Then, execute ``make OPT=opt LPS=none ZIMPL=false PARASCIP=true -j" to compile our software.

## Usage
You can solve _housing.linereg in ``data" by the following command

 - ./bin/scip -f data/housing.linereg -s settings/default.set

## Definition of Data
Our software reads from the second line of data file.
The definition of our data is as follows:
- [2nd line] the number of data
- [3rd line] the number of the explanatory variables
- [4th line] the index of the response variable(Begining is 1)
- [from 5th line] your data

See data files in the directory data for more details.

