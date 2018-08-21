#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   mode=$2

   echo "../bin/GENERATE_LPFILE $file $mode"
   ../bin/GENERATE_LPFILE $file $mode
   mv cplex.lp ${name}.lp
}

dirname=$1
mode=$2

for file in `ls ${dirname}/*.dat`
do
   run $file $time $mode
done

