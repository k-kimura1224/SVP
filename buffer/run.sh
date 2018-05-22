#/bin/sh
#chmod u+x gen.sh

thread=$1
time=$2

cd ../Split_miqp/run
./clean.sh
./run.sh ../../problems/TEST $time 2

cd ../../LPFILES/BKZ20LP_TEST
./clean.sh
./run.sh $time

cd ../../SVPSOLVER/run
./clean.sh
./run.sh ../../problems/TEST $thread $time

cd ../../buffer

echo "Done!!"
