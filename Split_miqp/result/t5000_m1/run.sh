#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   time=$2
   mode=$3

   echo "../bin/GENERATE_SPLIT_MIQP $file $time $mode > ${name}.log"
   ../bin/GENERATE_SPLIT_MIQP $file $time $mode > ${name}.log
}

dirname=$1
time=$2
mode=$3

for file in `ls ${dirname}/*.dat`
do
   run $file $time $mode
done

