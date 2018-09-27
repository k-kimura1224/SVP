#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   thread=$2
   time=$3

	export OMP_NUM_THREADS=1
   echo "../../bin/SVPSOLVER.opt -f $file -p $thread -t $time -q > ${name}.log"
   ../../bin/SVPSOLVER.opt -f $file -p $thread -t $time -q > ${name}.log
}

dirname=$1
thread=$2
time=$3

for file in `ls ${dirname}/*.dat`
do
   run $file $thread $3
done

