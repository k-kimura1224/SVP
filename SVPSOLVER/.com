alias opt42="./bin/SVPSOLVER.opt -f ../problems/TEST/BKZ20-SVP_n42_seed47.dat -p 1 -m 8 -e"
alias dbg42="./bin/SVPSOLVER.dbg -f ../problems/TEST/BKZ20-SVP_n42_seed47.dat -p 1 -m 8 -e"
