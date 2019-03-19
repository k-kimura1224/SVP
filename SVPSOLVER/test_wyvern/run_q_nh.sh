#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   opt=$2
   thread=$3
   time=$4
   mem=$5

	export OMP_NUM_THREADS=1
   echo "../bin/SVPSOLVER.${opt} -f $file -p $thread -t $time -m $mem -q -H -1 > ${name}.log"
   ../bin/SVPSOLVER.${opt} -f $file -p $thread -t $time -m $mem -q -H -1 > ${name}.log
}

dirname=$1
opt=$2
thread=$3
time=$4
mem=$5

for file in `ls ${dirname}/*.dat`
do
   run $file $opt $thread $time $mem
done

