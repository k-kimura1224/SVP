#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   thread=4
   time=86400

	export OMP_NUM_THREADS=1
   echo "../bin/SVPSOLVER.dbg -f $file -p $thread -t $time > ${name}.log"
   ../bin/SVPSOLVER.dbg -f $file -p $thread -t $time > ${name}.log
}

dirname=$1

for file in `ls ${dirname}/*.dat`
do
   run $file
done

