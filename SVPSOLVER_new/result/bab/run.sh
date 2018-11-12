#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   thread=$2
   time=$3
	memory=$4

	export OMP_NUM_THREADS=1
   echo "../../bin/SVPSOLVER.opt -f $file -p $thread -t $time -q -m $memory -r > ${name}.log"
   ../../bin/SVPSOLVER.opt -f $file -p $thread -t $time -q -m $memory -r > ${name}.log
}

dirname=$1
thread=$2
time=$3

#for file in `ls ${dirname}/*.dat`
#do
#   run $file $thread $3
#done
#for file in `ls ../../../link/*/*.dat`
#do
#   run $file 32 86400
#done

#for i in 40 43 46 49 52 55 58 61
for i in 52 55
do
	for j in 0 1 2 3 4
	do
		file=../../../link/B20_${i}s/BKZ20-SVP_n${i}_seed${j}.dat
   	run $file 32 86400 128
	done
done




